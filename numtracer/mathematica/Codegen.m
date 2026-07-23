(* ::Package:: *)
(* Code generation. NumTracer owns ONLY the tensor part — contracting each
   component numerically to a polynomial (MPoly) that is Horner-lowered to a real
   straight-line kernel. Everything scalar (the dressing/regulator coefficients, CSE,
   powr<n>, the function/class/header boilerplate, clang-format, write-if-changed)
   is delegated to FunKit's mature COEN emitter (CppForm / MakeCppFunction /
   MakeCppClass / MakeCppHeader / WriteCodeToFile), which produces the flat
   straight-line kernel form.

   The seam: each component's scalar result is bound to a C++ identifier in the
   kernel-body PREAMBLE we emit; that identifier appears as a *string placeholder*
   in the Mathematica integrand `Σ coeff_i × trace_i`, which FunKit's CppForm emits
   verbatim (CExpression[a_String] := a) while CSE-ing the coefficients around it. *)
(* Generation goes through the numeric matrix-product backend (`builderInv`/`compileTInv` + the
   generator below); generated kernels are validated against FormTracer (FORM) oracles in the
   test suite. *)
(* ============================================================================ *)
(* ==== INVARIANT-BASIS path: contract in scalar-product symbols, not frame  === *)
(* ==== components. The Lorentz trace of each diagram is emitted as an `et::inv` *)
(* ==== network into a build-time GENERATOR (gen_<name>_inv.cpp); the generator *)
(* ==== runs reduce->cancel->lower at codegen time (heap reclaimed, ~FORM op    *)
(* ==== count) and prints a committed straight-line `<name>_inv_kernels.hh`.    *)
(* ==== The kernel just fills the few fundamental symbols and calls trN(f).     *)
(* ==== This closes the A3 runtime gap that the frame-component basis cannot.   *)
(* ============================================================================ *)
(* Verbose-diagnostics gate. The profiling / CSE / probe / timing traces below ([prof], [cse],
   [probe], [diagpoly], [time]) are emitted through ntLog and stay SILENT unless this flag is set —
   so a normal generation run is quiet. Genuine "[NumTracer] ERROR" aborts and the "wrote:"/
   "unchanged:" file messages are always printed (plain Print). To see the diagnostics, set
   NumTracer`Private`$NumTracerVerbose = True before generating. ntLog evaluates its arguments only
   when verbose, so it must NEVER wrap a side-effecting computation — keep such work (e.g. an
   AbsoluteTiming around an assignment) in a With-binding and pass only the timing into ntLog. *)
(* Default silent, but env-controllable so a headless/CI run can turn the [cse]/[prof]/[time]
   diagnostics on without editing a .wls — mirrors the C++ side's NT_GEN_PROFILE. The density guard in
   tests/gen/regen_check.sh needs the [cse] sub-terms line, which is emitted through ntLog. *)
$NumTracerVerbose = (Environment["NT_GEN_VERBOSE"] =!= $Failed);

ntLog[args___] := If[TrueQ[$NumTracerVerbose],
    Print[args]];

(* WolframKernel resident set (MB) from /proc — the RSS the OOM killer sees, not MemoryInUse[]
   (allocated). Returns 0 where /proc is absent. *)
ntWolframRssMB[] := Quiet @ Check[
    Module[{s = ReadString["/proc/self/statm"]},
      If[StringQ[s], ToExpression[StringSplit[s][[2]]] * 4096. / 1048576., 0.]],
    0.];

(* a C++ double literal at full precision (exact rationals -> doubles). *)

cppNum[x_] := ToString[CForm[N[x, 17]]];

(* CppForm, flattened to a single line — the fill formulas are emitted INSIDE C++ string
   literals, so CppForm's line-wrapping (newlines) must be collapsed or they break the string. *)

cppFlat[e_] := StringReplace[FunKit`CppForm[e], {"\n" -> " ", "\r" -> " ", "\t" -> " "}];

(* ---- the single emission chokepoint --------------------------------------------------------
   EVERY generated file goes through here. Nothing else may call Export on generated source.

   The failure this exists for: an expression that reached the emitter WITHOUT being lowered to C++
   gets ToString'd verbatim into the file. It has happened three times — the ZAAqbq metric leak
   (`colFacG[ntMetric[...], <|...|>]`), a degenerate Gram emitting `return Indeterminate;`, and the
   four-quark Fierz flavour sum (`colFacG[Plus[...], <|...|>]`). Each time the symptom appeared
   layers away from the cause: a clang syntax error naming a column of a 7000-character line, or —
   the dangerous variant — a leak whose text happens to BE valid C++ and compiles into a silently
   wrong kernel.

   A textual assertion catches the whole class at the one place it must pass through, regardless of
   which upstream dispatcher grew a hole. Structural guards upstream (tleak/colleak/eagernn) give
   better messages and should stay; this is the backstop that cannot be forgotten.

   `nt` heads are matched with a word boundary so ordinary C++ identifiers containing "nt"
   (`constant`, `int`, `point`) do not trip it. *)

$ntCppLeakPatterns =
  {
    RegularExpression["(?<![A-Za-z0-9_])nt[A-Z][A-Za-z0-9]*\\["],
    (* any nt* DSL head, un-lowered *)
    "colFacG[",
    "compileColG[",
    "compileTInv[",
    "colourFactorProd[",
    "Indeterminate",
    "DirectedInfinity",
    "ComplexInfinity",
    "$Failed",
    "Missing["};

ntExportCpp[file_, text_] := (
    Function[patt,
        Module[{pos = StringPosition[text, patt, 1]},
          If[pos =!= {},
(* quote the leak IN CONTEXT — the matched token alone rarely identifies which structure it
   came from, and the emitted lines are thousands of characters wide. *)
            Message[MakeNTKernel::cppleak, "ntExportCpp", StringTake[text, First[pos]], file <> "\n  context: ..." <> StringTake[text, {Max[1, pos[[1, 1]] - 150], Min[StringLength[text], pos[[1, 2]] + 150]}] <> "..."];
            Abort[]]]
      ] /@ $ntCppLeakPatterns;
    Export[file, text, "Text"]);

(* ---- GlobalCollect: dressing-coefficient decomposition (the Route-B front-end) ----------------
   A diagram's dressing/kinematic coefficient is decomposed into a SUM of monomials, each
   {numericCoeff, {atomExpr...}}. An ATOM is a maximal non-numeric multiplicative factor that the
   kernel evaluates ONCE (a propagator dressing `ZA3[...]`, a regulator `RB[...]`, a composite
   denominator `Power[D,-1]`, a kinematic factor `cos1`/`Sqrt[...]`). Power[b,n] with a positive
   integer n expands to n copies of atom b (so it lowers as b^n); any other factor — including a
   negative/fractional power or a function call — is one atom. The numeric part folds into the
   collector coefficient; the atoms become inert DRESSING SYMBOLS (a `dr` tag) so identical
   (scalar-product × dressing) monomials from different diagrams MERGE — the cross-diagram collection
   FORM does. Returns {monomials, ...} where each monomial is {num (machine real), {atomExpr...}}. *)

dressMonomials[expr_] := Module[{terms},
    terms = Expand[expr];
    terms =
      If[Head[terms] === Plus,
        List @@ terms,
        {terms}];
    Function[term,
        Module[{num = 1, atoms = {}, factors},
          factors =
            If[Head[term] === Times,
              List @@ term,
              {term}];
          Do[
            Which[
              NumberQ[f],
                num *= f,
              MatchQ[f, Power[_, _Integer?Positive]],
                atoms = Join[atoms, ConstantArray[First[f], Last[f]]],
              True,
                AppendTo[atoms, f]],
            {f, factors}];
          {N[num, 17], atoms}]
      ] /@ terms];

(* ---- inv dialect: same tensor structure as `compileT`, but emitting the
        et::inv builder calls the generator's wrapper templates expose
        (tproj / lmetric / lvec / sc / contract / add). *)

builderInv[ntMetric[mu_, nu_], ids_, env_, mask_, nc_] :=
  "lmetric<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ">()";

builderInv[ntVec[q_, mu_], ids_, env_, mask_, nc_] :=
  "lvec<" <> ToString[ids[mu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ">()";

builderInv[ntTransProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "tproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ">()";

builderInv[ntLongProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "lproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ">()";

builderInv[ntMagneticProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "mproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["InvS"]] <> ">()";

builderInv[ntElectricProj[q_, mu_, nu_], ids_, env_, mask_, nc_] :=
  "eproj<" <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[mask[q]] <> ", " <> ToString[env[q]["Inv"]] <> ", " <> ToString[env[q]["InvS"]] <> ">()";

builderInv[ntEpsilon[a_, b_, c_, d_], ids_, env_, mask_, nc_] :=
  "leps<" <> ToString[ids[a]] <> ", " <> ToString[ids[b]] <> ", " <> ToString[ids[c]] <> ", " <> ToString[ids[d]] <> ">()";

(* A Lorentz factor as a single `network::Elem{...}` literal (for a collected Dirac slot's per-option
   `netFacs`, Stage 4). Mirrors builderInv's id/momentum/atom resolution but emits the Elem aggregate the
   numeric backend appends to the net, rather than a NetVal builder. Field order (network.hpp):
   {kind, a, b, vid, inv, vlc, c, d, invS}. A projector's momentum rides `vid = env Base` (elem_to_nelem
   reconstructs it as {{1.0, vid}}); a vector's rides `vlc`. *)
builderInvElem[ntMetric[mu_, nu_], ids_, env_] :=
  "Elem{Elem::Metric, " <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", -1, -1, {}}";
builderInvElem[ntVec[q_, mu_], ids_, env_] :=
  "Elem{Elem::Vector, " <> ToString[ids[mu]] <> ", -1, -1, -1, {{1.0, " <> ToString[env[q]["Base"]] <> "}}}";
builderInvElem[ntTransProj[q_, mu_, nu_], ids_, env_] :=
  "Elem{Elem::ProjT, " <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[env[q]["Inv"]] <> ", {}}";
builderInvElem[ntLongProj[q_, mu_, nu_], ids_, env_] :=
  "Elem{Elem::ProjL, " <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[env[q]["Inv"]] <> ", {}}";
builderInvElem[ntMagneticProj[q_, mu_, nu_], ids_, env_] :=
  "Elem{Elem::ProjM, " <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", -1, {}, 0, 0, " <> ToString[env[q]["InvS"]] <> "}";
builderInvElem[ntElectricProj[q_, mu_, nu_], ids_, env_] :=
  "Elem{Elem::ProjE, " <> ToString[ids[mu]] <> ", " <> ToString[ids[nu]] <> ", " <> ToString[env[q]["Base"]] <> ", " <> ToString[env[q]["Inv"]] <> ", {}, 0, 0, " <> ToString[env[q]["InvS"]] <> "}";
builderInvElem[ntEpsilon[a_, b_, c_, d_], ids_, env_] :=
  "Elem{Elem::Epsilon, " <> ToString[ids[a]] <> ", " <> ToString[ids[b]] <> ", -1, -1, {}, " <> ToString[ids[c]] <> ", " <> ToString[ids[d]] <> "}";

scaleStrInv[str_, 1] := str;

scaleStrInv[str_, 1.] := str;

scaleStrInv[str_, s_] := "sc<numtracer::Lit<numtracer::Cx{" <> cppNum[s] <> ", 0.0}>>(" <> str <> ")";

wrapContractInv[{one_}] := one;

wrapContractInv[many_] := "contract(" <> StringRiffle[many, ", "] <> ")";

MakeNTKernel::eagernn = "compileTInv: an eagerly-summed structure has a NON-NUMERIC per-structure scalar coefficient, which et::add cannot fold (it scales each structure by a compile-time Lit). The sum should have been distributed by expandBridges (DSL.m) or collected into an ntDressedNum. Offending sum:\n`1`";

MakeNTKernel::tleak = "compileTInv: un-lowered TENSOR structure reached the scalar fallthrough — it would be CForm'd into C++ as a bare scalar with its indices silently dropped (this is what caused the ZAAqbq metric leak). Every tensor head must be handled by builderInv or one of the Power/Times/Plus branches. Offending structure:\n`1`";

MakeNTKernel::colleak = "compileColG: a factor of a CONSTANT SU(N) component matched none of the colFacG head rules, so it would be ToString'd into the generated C++ as raw Mathematica (e.g. `colFacG[Plus[...], <|...|>]`), which does not compile. Every factor of a constant component must be one of the six group heads (ntSUNf/ntSUNDeltaAdj/ntSUNT/ntSUNDeltaFund/ntSUNDiag{Fund,Adj}). A Plus here means a colour/flavour sum reached the flat compiler instead of compileColGSum; any other head means the \"Constant\" head list in DSL.m analyseDiagram has drifted. Offending factor:\n`1`";

(* NB: no backquoted code fragments in this string — a backquoted word is a StringForm SLOT, so
   quoting an identifier that way makes the message itself fail to format (StringForm::sfr). *)

MakeNTKernel::colrest = "splitColourGroupsInv: after splitting a branch into its colour product and its Lorentz/Dirac remainder, SU(N) head(s) are still present in the REMAINDER. The remainder goes to the Lorentz-only builderInv/compileDirac, which has no rule for a group head, so it would be CForm'd into the generated C++ as raw Mathematica (builderInv[ntSUNDeltaFund[...], <|...|>]). The split is Cases/DeleteCases at LEVEL 1, so it only sees group heads that are BARE factors — one buried inside a Power (a CLOSED colour loop: deltaFund[N,i,j]^2) or a surviving Plus slips through. Offending head(s):\n`1`\nRemainder:\n`2`";

MakeNTKernel::colpow = "compileColG: a colour/flavour SUM raised to the integer power `1`. Expanding it by repetition would duplicate the summands' index labels, so the same label would appear on 2k tensors and the et/SUNNet contraction would silently mis-pair them into a wrong number — and checkLabels has already run by this point, so nothing downstream would notice. (This is the colour analogue of NumTrace::bridgepow.) Refusing instead. Offending base:\n`2`";

MakeNTKernel::cppleak = "`1`: the generated source still contains un-lowered Mathematica — the text `2` appears in it. Writing it would produce a file that either fails to compile or, worse, compiles into a silently wrong kernel. This means some expression reached the emitter without being turned into C++; the fragment above should identify which. Offending file:\n`3`";

(* MEMOIZED: the projector/Lorentz net builder is called with mostly repeated inputs — the same
   transverse-projector structures recur across every diagram/branch, so a dense flow makes orders of
   magnitude more calls than it has distinct arguments. Cache by a hash of the args; $ctCache is cleared
   per generation in mkGenerateKernel. Output-preserving (a pure function of its inputs); the recursion
   is memoised too. *)

compileTInv[e_, ids_, env_, mask_, nc_] := With[{h = Hash[{e, ids, env, mask, nc}]},
    Lookup[$ctCache, h, $ctCache[h] = compileTInvBody[e, ids, env, mask, nc]]];

compileTInvBody[e_, ids_, env_, mask_, nc_] := Which[
    tensorQ[e],
      {builderInv[e, ids, env, mask, nc], 1},
    Head[e] === Power && IntegerQ[e[[2]]] && e[[2]] >= 1 && !scalarQ[e],
(* A TENSOR raised to an integer power = that many copies sharing the SAME index labels, i.e. a
   closed self-contraction (e.g. ntMetric[v1,v2]^2 = g_{v1 v2} g_{v1 v2} = D). Expand into n
   contracted copies so the numeric index-elimination folds it to a number. Without this it falls
   through to the scalar branch below and is CForm'd into undeclared C++ (the ZAAqbq metric leak). *)Module[{cs = Table[compileTInv[e[[1]], ids, env, mask, nc], {e[[2]]}]},
        {wrapContractInv[cs[[All, 1]]], Times @@ cs[[All, 2]]}],
    Head[e] === Times,
      Module[{parts = List @@ e, sc, tn, cs},
        sc = Select[parts, scalarQ];
        tn = Select[parts, !scalarQ[#]&];
        cs = compileTInv[#, ids, env, mask, nc]& /@ orderFactors[tn];
        {wrapContractInv[cs[[All, 1]]], (Times @@ sc) (Times @@ cs[[All, 2]])}],
    Head[e] === Plus,
      Module[{cs = compileTInv[#, ids, env, mask, nc]& /@ (List @@ e)},
        If[!AllTrue[cs[[All, 2]], NumericQ],
          Message[MakeNTKernel::eagernn, e];
          Abort[]];
        {"add(" <> StringRiffle[MapThread[scaleStrInv, {cs[[All, 1]], cs[[All, 2]]}], ", "] <> ")", 1}],
    (* a genuine scalar coefficient: no builder, the expression IS the scalar *)scalarQ[e],
      {"", e},
(* Anything else still carries a TENSOR head but matched none of the branches above, so it
   would be emitted as a bare C++ scalar with its indices silently dropped (the ZAAqbq metric
   leak, see the Power comment above). Fail loudly instead — this is the Lorentz/colour
   analogue of the Dirac leak guard below. *)True,
      Message[MakeNTKernel::tleak, Short[e, 6]];
      Abort[]];

(* ---- colour/Lorentz sector split for a colour-ENTANGLED Lorentz/Dirac component -----------
   The full quark-gluon vertex basis (AqbqDirect147 structures 4/7) produces a quark-loop
   component that is a SUM whose terms pair different colour orderings (T^{c1}T^{c2} vs
   T^{c2}T^{c1}) with different Dirac traces — colour and Dirac do NOT factor globally, only per
   term. expandDiracComponent then leaves the colour tensors inside the (would-be Lorentz) result,
   where the Lorentz-only builderInv cannot lower them. The fix: group the Dirac-traced expression
   by its colour-factor product; each group has ONE colour structure (folded numerically as a
   SUNNet — exactly like a constant adjoint/fundamental component) times a pure-Lorentz
   polynomial (the inv net). The diagram's trace is Sum_group colv_group * lorentzPoly_group,
   emitted as several (sunNet, lorentzNet) entries that the per-diagram combination already sums. *)

ctHeadsInv = {_ntSUNf, _ntSUNDeltaAdj, _ntSUNT, _ntSUNDeltaFund, _ntSUNDiagFund, _ntSUNDiagAdj};

colourEntangledQ[e_] := !FreeQ[e, Alternatives @@ ctHeadsInv];

(* A group head OR an integer power of one. Splitting a term into "colour" and "the Lorentz/Dirac
   rest" is a LEVEL-1 Cases/DeleteCases over the factor list, so it matches only what is a BARE
   factor — and a CLOSED colour/flavour loop collapses two identical deltas into deltaFund[N,i,j]^2
   (= N), which is a Power, not a head. Matched by the head pattern alone, such a factor is silently
   left in the remainder and handed to the Lorentz-only builderInv, which CForm's it into the
   generated C++. compileColG already knows how to expand this power into repeated SUNNet factors;
   it just has to be COLLECTED here first. (The four-quark Fierz gate reached exactly this: the
   epsilon-pair expansion produces delta products, and a closed flavour loop squares one of them.) *)

ctFacInv = Alternatives[Alternatives @@ ctHeadsInv, Power[Alternatives @@ ctHeadsInv, _Integer?Positive]];

colourFactorProd[term_] := Times @@
    Cases[
      If[Head[term] === Times,
        List @@ term,
        {term}],
      ctFacInv];

mergeColNet["SUNNet{}", b_] := b;

mergeColNet[a_, "SUNNet{}"] := a;

mergeColNet[a_, b_] :=
  "SUNNet{" <> StringDrop[StringDrop[a, 7], -1] <> "," <> StringDrop[StringDrop[b, 7], -1] <> "}";

(* ---- per-component diagonal-dressing registry (ntSUNDiag{Fund,Adj}) ------------------------
   A diag-dressed group δ dresses SELECTED components with distinctly-named scalar dressings (e.g.
   the Cartan directions of a condensate) and DROPS the rest. Its `spec` is a rules list
   {c1 -> name1, …, Default -> defName}: `ci` are 1-based component indices, `namei` scalar dressing
   symbols; components with no rule (and no Default) vanish. colFacG registers each distinct
   (name, scale) leaf under a small integer id `dr` and bakes a per-component id vector
   (component → dr, -1 = drop) into the emitted sun<n>.diag{Fund,Adj}(...,{d0,…}) factor. The C++ seam
   folds the net to a SUNPoly over these ids, and the integrand multiplies in the runtime sum
   Σ_t coeff_t Π name(scale) — an ordinary scalar-dressing token, no array. Reset per generation. *)

$diagDrTable = <||>; $diagDrByKey = <||>; $diagDrCounter = 0;

resetDiagDr[] := (
    $diagDrTable = <||>;
    $diagDrByKey = <||>;
    $diagDrCounter = 0;);

diagDrId[name_, scale_] := Module[{key = {name, scale}},
    If[!KeyExistsQ[$diagDrByKey, key],
      $diagDrByKey[key] = $diagDrCounter;
      $diagDrTable[$diagDrCounter] = <|"Name" -> name, "Scale" -> scale|>;
      $diagDrCounter++];
    $diagDrByKey[key]];

(* Parse a diag-dressing spec into a per-component id vector of length `dim` (component 0..dim-1;
   1-based physics index = v+1). Named components get diagDrId[name, scale]; a Default -> name rule
   fills the rest; unmatched components are -1 (dropped). *)

diagComp2Dr[spec_, scale_, dim_] := Module[{rules = Flatten[{spec}], named, def},
    named = Association[Cases[rules, (c_Integer -> nm_) :> (c -> diagDrId[nm, scale])]];
    def = Cases[rules, (Default -> nm_) :> diagDrId[nm, scale]];
    def =
      If[def === {},
        -1,
        First[def]];
    Table[Lookup[named, v + 1, def], {v, 0, dim - 1}]];

diagVecStr[vec_] := "{" <> StringRiffle[ToString /@ vec, ","] <> "}";

(* ---- scalar-dressing registry (symbolic dressing collection: ntDressedNum slots) -------------
   Each DISTINCT dressing atom — a maximal non-numeric multiplicative factor in a dressed numerator's
   per-structure coefficient (a propagator dressing `Zq[...]`, `hSigL[...]`, a regulator `RB[...]`, a
   composite denominator `Power[D,-1]`, a kinematic `cos1`) — is interned under a small integer id baked
   into the emitted DSlotOpt and into the generator's fm.dress table. So identical atoms across slots /
   diagrams share ONE `f[]` slot and the runtime evaluates each dressing call ONCE (the cross-diagram
   collection FORM does). The atom expression is already frame-resolved (ntSP/ntVec components
   substituted) so its C++ fill is `cppFlat[atom]`. Reset per generation alongside resetDiagDr. *)

$drTable = <||>; $drByKey = <||>; $drCounter = 0;

resetDr[] := (
    $drTable = <||>;
    $drByKey = <||>;
    $drCounter = 0;);

drAtomId[atom_] := Module[{key = atom},
    If[!KeyExistsQ[$drByKey, key],
      $drByKey[key] = $drCounter;
      $drTable[$drCounter] = atom;
      $drCounter++];
    $drByKey[key]];

(* Decompose a (frame-resolved) dressed-structure coefficient into {Cx numeric, {drAtomId…}}: numbers
   fold into the complex coefficient; a positive-integer power b^n expands to n atom copies; every other
   non-numeric factor is one atom (interned via drAtomId). Mirrors dressMonomials but interns atoms. *)

drDecompose[coeff_] := Module[{
    factors =
      If[Head[coeff] === Times,
        List @@ coeff,
        {coeff}],
    num = 1,
    ids = {}},
    Do[
      Which[
        NumberQ[f],
          num *= f,
        MatchQ[f, Power[_, _Integer?Positive]],
          ids = Join[ids, ConstantArray[drAtomId[First[f]], Last[f]]],
        True,
          AppendTo[ids, drAtomId[f]]],
      {f, factors}];
    {N[num, 17], Sort[ids]}];

(* Chunk a Lorentz polynomial into several inv nets of <= $ntInvChunk top-level terms each, so no
   single generated net-builder function becomes a giant nested add() that blows up the g++ -O0
   compile (a quark box with the full quark-gluon vertex basis can be tens of thousands of nodes in
   ONE net). Returns a list of {lorentzNetString, scalar}; the per-diagram combination sums the
   chunks (same colour + coeff) into one trace, so chunking is transparent. Zero -> dropped; a pure
   number -> a constant net (konst). *)

$ntInvChunk = 150;

(* max chars of net-builder elements packed into ONE emitted function (see ntChunkDef below). A net
   builder is normally emitted as a single braced-init list of its sub-term elements:
     std::vector<std::vector<DChainTok>> dch<i>(){ return {E1, ..., En}; }
   On a dense flow that one expression grows to megabytes / tens of thousands of elements. At -O0 the
   back end handles it as ONE basic block with thousands of live temporaries, so a single such TU can
   cost minutes and several GB of RSS — the direct cause of both a very long compile AND an OOM, since
   several such units run in parallel. Worse, the braced-init materialises every element temporary in
   one full-expression, so a big enough builder overflows the default stack at RUNTIME. Chunking the
   element list into size-bounded helper functions fixes all three, and — since each helper is its own
   top-level def — lets the unit bin-packer spread them across TUs, so no single unit carries a whole
   giant net. Kept well under $ntUnitChars so chunks stay packable. *)

$ntDefChunk = 60000;

(* target chars per generator TU, and the hard cap on how many TUs a flow may split into. The cap must
   stay above what the size target asks for, or units silently grow back past the target. *)

$ntUnitChars = 250000;

$ntUnitCap = 512;

(* per-job RAM estimate (GB) used to bound the parallel compile; NT_GEN_JOB_MEM overrides. With defs
   chunked, a unit's peak RSS is a few hundred MB, so 1.5 is deliberately conservative headroom. *)

$ntGenJobMemGB =
  With[{e = Environment["NT_GEN_JOB_MEM"]},
    If[StringQ[e] && NumericQ[Quiet @ ToExpression[e]] && ToExpression[e] > 0,
      N @ ToExpression[e],
      1.5]];

(* MemAvailable: what the kernel can hand out without swapping. Read fresh (not at package load) —
   the Wolfram kernel that just built the flow may itself be holding many GB. *)

ntAvailMemGB[] := Quiet @
    Check[
      Module[{m = StringCases[Import["/proc/meminfo", "Text"], "MemAvailable:" ~~ Whitespace ~~ d : DigitCharacter.. ~~ Whitespace ~~ "kB" :> ToExpression[d], 1]},
        If[m === {},
          N[MemoryAvailable[] / 2^30],
          First[m] / 1024. / 1024.]],
      $Failed];

(* max concurrent cxx processes in the generator compile. Peak RAM ~ jobs x per-unit RSS, so this is
   bounded by BOTH cores and free memory: a cores-only cap will OOM the box on a dense flow, where a
   single unit can take GBs. Def chunking now bounds a unit's size, so the memory term is belt-and-braces
   that keeps a future monster flow from thrashing instead of failing loudly. NT_GEN_JOBS pins the count
   outright; NT_GEN_JOB_MEM tunes the per-job estimate. Delayed (:=) so the memory reading is taken at
   compile time, not at package load. *)

ntCompileJobs[] := With[{e = Environment["NT_GEN_JOBS"], avail = ntAvailMemGB[]},
    Which[
      StringQ[e] && IntegerQ[Quiet @ ToExpression[e]] && ToExpression[e] > 0,
        ToExpression[e],
      IntegerQ[$ProcessorCount],
        Max[
          2,
          Min[
            24,
            $ProcessorCount,
            If[NumericQ[avail] && avail > 0,
              Floor[avail / $ntGenJobMemGB],
              24]]],
      True,
        4]];

$ntCompileJobs := ntCompileJobs[];

(* Emit ONE net builder, chunking it when its element list is too big to sit in a single function.
   Small builders keep the old single-def form byte-for-byte. Big ones become

     void dch<i>_c0(std::vector<std::vector<DChainTok>>& o){ o.push_back(E1); o.push_back(E2); ... }
     ...
     std::vector<std::vector<DChainTok>> dch<i>(){
       std::vector<std::vector<DChainTok>> o; o.reserve(n); dch<i>_c0(o); ...; return o; }

   Order is preserved, so the assembled vector is element-for-element identical to the braced-init
   one — and strictly less work at runtime, since a braced-init-list copies every element while
   push_back moves the temporary. Returns {defs, decls}: the helpers go into `allDefs` as ordinary
   top-level defs (so the bin-packer may scatter them across units), and their forward declarations
   go into the shared decl header, which is where the cross-unit calls resolve. *)

ntChunkDef[name_String, ret_String, elems_List] :=
  If[Length[elems] <= 1 || Total[StringLength /@ elems] <= $ntDefChunk,
    {{ret <> " " <> name <> "(){ return {" <> StringRiffle[elems, ", "] <> "}; }"}, ""},
    Module[
      {cs, chunks, nc, defs},
      (* cumulative chars / chunk size is nondecreasing, so equal keys form contiguous runs *)
      cs = Ceiling[Accumulate[(StringLength /@ elems) + 2] / $ntDefChunk];
      chunks = SplitBy[Transpose[{elems, cs}], Last][[All, All, 1]];
      nc = Length[chunks];
      defs = MapIndexed["void " <> name <> "_c" <> ToString[#2[[1]] - 1] <> "(" <> ret <> "& o){ " <> StringJoin[("o.push_back(" <> # <> "); ")& /@ #1] <> "}"&, chunks];
      {Append[defs, ret <> " " <> name <> "(){ " <> ret <> " o; o.reserve(" <> ToString[Length[elems]] <> "); " <> StringJoin[Table[name <> "_c" <> ToString[k - 1] <> "(o); ", {k, nc}]] <> "return o; }"], StringJoin[Table["void " <> name <> "_c" <> ToString[k - 1] <> "(" <> ret <> "&);\n", {k, nc}]]}
    ]];

(* the same over every net of one family (dnet/lnet/dch/dsl); returns {flatDefs, declString}. *)

ntChunkDefs[prefix_String, ret_String, elemLists_List] := If[elemLists === {},
    {{}, ""},
    Module[{r = MapIndexed[ntChunkDef[prefix <> ToString[#2[[1]] - 1], ret, #1]&, elemLists]},
      {Flatten[r[[All, 1]]], StringJoin[r[[All, 2]]]}]];

(* C++ compiler for the build-time generator. The emitted generator is ordinary numeric C++, so
   either g++ or clang++ compiles it with the same flags (-std=c++20 -ftemplate-depth=4000
   -O2/-O0 -fno-exceptions -fno-rtti -pthread -I -c -o). Resolution order: the NT_GEN_CXX env
   override (verbatim), else PREFER clang++ when present (it compiles the -O0 net-builder units
   ~2.5x faster than g++ — u33 5.8 s -> 2.3 s), else a compiler detected via CCompilerDriver
   (GCC -> g++, Clang -> clang++; binary taken from the driver's CompilerInstallation dir), else
   "g++" on PATH. clang objects link fine against the g++-built libNumTracer.a (shared libstdc++
   ABI on Linux). *)

resolveGenCxx[] := Module[{env = Environment["NT_GEN_CXX"], path, comps, clangPick, pick, name, inst, exe, full},
    If[StringQ[env] && StringTrim[env] =!= "",
      Return[StringTrim[env]]];
    Quiet @ Needs["CCompilerDriver`"];
    comps = Quiet @ Check[CCompilers[], {}];
    If[!ListQ[comps],
      comps = {}];
    comps = Select[comps, AssociationQ[#] && StringQ[Lookup[#, "CompilerInstallation"]]&];
    (* 1) prefer a Clang known to CCompilerDriver *)
    clangPick = SelectFirst[comps, StringContainsQ[ToString @ Lookup[#, "Name", ""], "Clang", IgnoreCase -> True]&, None];
    If[clangPick =!= None,
      full = FileNameJoin[{Lookup[clangPick, "CompilerInstallation"], "clang++"}];
      Return[
        If[FileExistsQ[full],
          full,
          "clang++"]]];
    (* 2) prefer clang++ on PATH even if CCompilerDriver did not enumerate it *)
    path = Environment["PATH"];
    If[StringQ[path] && AnyTrue[StringSplit[path, ":"], FileExistsQ[FileNameJoin[{#, "clang++"}]]&],
      Return["clang++"]];
    (* 3) fall back to the CCompilerDriver default (typically g++) *)
    If[comps === {},
      Return["g++"]];
    pick = SelectFirst[comps, Lookup[#, "Compiler"] === Quiet @ DefaultCCompiler[]&, First[comps]];
    name = ToString @ Lookup[pick, "Name", ""];
    inst = Lookup[pick, "CompilerInstallation"];
    exe =
      If[StringContainsQ[name, "Clang", IgnoreCase -> True],
        "clang++",
        "g++"];
    full = FileNameJoin[{inst, exe}];
    If[FileExistsQ[full],
      full,
      exe]];

chunkLorInv[lorExpr_, ids_, env_, mask_, nc_] := Which[
    lorExpr === 0,
      {},
    scalarQ[lorExpr],
      {{"konst(" <> cppNum[lorExpr] <> ")", 1}},
    True,
      Module[{
        terms =
          If[Head[lorExpr] === Plus,
            List @@ lorExpr,
            {lorExpr}]},
        (compileTInv[Total[#], ids, env, mask, nc]& /@ Partition[terms, UpTo[$ntInvChunk]])]];

(* Build inv nets for keepLor * inner, where keepLor is the big COMMON Lorentz polynomial (the
   projector x kinematics, shared by every branch of a colour group) and inner is the per-group
   Dirac-trace sum. keepLor's terms are chunked (<= $ntInvChunk each) and inner is multiplied onto each
   chunk, so the projector net is built ONCE PER CHUNK instead of once per branch. Building the
   projector net per branch (the old `chunkLorInv` over a fully-multiplied-out lorSum) was the
   string-build hot spot — a single 3-/4-point quark component could sit minutes inside compileTInv. *)

chunkLorProd[keepLor_, inner_, ids_, env_, mask_, nc_] := Which[
    inner === 0 || inner === 0.,
      {},
    keepLor === 1,
      chunkLorInv[inner, ids, env, mask, nc],
    Head[keepLor] === Plus,
      (compileTInv[Total[#] inner, ids, env, mask, nc]& /@ Partition[List @@ keepLor, UpTo[$ntInvChunk]]),
    True,
      {compileTInv[keepLor inner, ids, env, mask, nc]}];

(* {colourNet, lorentzNet, scalar} entries per colour-structure group of a colour-entangled component.
   DISTRIBUTE only the entangled Pluses (the vertex colour-ordering / commutator sub-sums, small) — NOT
   the big pure-Lorentz angular polynomial, which stays factored in each branch. Each branch is split
   into its colour-factor product (folded numerically → SUNNet) and its colour-free rest; the rest
   is emitted as a `dirac_value` net (gammas traced in C++) when it carries a gamma chain, else as a
   chunked Lorentz net. Branches are grouped by colour structure and summed. The Dirac trace is thus
   contracted in the C++ generator, not expanded symbolically here. *)
(* ---- struct-7 σ^{μν} folding: keep the bare γ-commutator as ONE token ------------------------
   The quark-gluon-vertex struct-7 tensor σ^{μν}=(i/2)[γ^μ,γ^ν] arrives from FunKit as a bare 2-term
   antisymmetric γ-pair `Plus`  s·(γ(X)γ(Y) − γ(Y)γ(X))  (FunKit has no σ primitive). Left alone,
   splitColourGroupsInv's `Expand` distributes it into TWO full Dirac traces. `foldDiracSigma`
   collapses that Plus into a single `ntSigma[legA, legB, din, dout]` token (each leg a slashed
   momentum {"slash",mom} or a free gluon id {"free",mu}), so the commutator is traced ONCE (the C++
   engine folds [A,B] as a block-diagonal 2×2 factor — `dcomm*` in et/inv/dirac.hpp). The i/2 and any
   sign live in the SCALAR (the Plus is already a bare bracket). It is a pure OPTIMIZATION: when the
   Plus is not an UNAMBIGUOUS commutator the recognizer returns $Failed and it distributes as before. *)
(* the (type,value) leg of a gamma `g` within a term's factor list `tf`. The gamma's Lorentz label μ
   is carried by the factor(s) holding ntVec[_,μ): either a single ntVec[mom,μ] OR a nested momentum
   sum like (ntVec[l1,μ]−ntVec[p2,μ]) — the LOOP σ legs are momentum linear combinations (l1−p2, …),
   not single vecs. The leg is returned as a SORTED list of {coeff, q} pairs, each q the LITERAL ntVec
   momentum (an env key by construction: momentumOf collected every ntVec momentum into the env basis),
   so legStr maps it straight to a vlc without re-decomposing over a non-atomic basis. A μ carried by
   no ntVec ⇒ a free (open) gluon id {"free",μ}; an ill-formed leg (a μ-bearing factor that is not a
   plain coeff·ntVec sum) ⇒ $Failed (stay conservative — the commutator then distributes as before). *)

ntSigmaLeg[leg_, termFactors_] := Module[{mu = First[leg], momFactors, terms},
    momFactors = Times @@ Select[termFactors, !FreeQ[#, ntVec[_, mu]]&];
    If[momFactors === 1,
      Return[{"free", mu}]];
    terms =
      If[Head[momFactors] === Plus,
        List @@ momFactors,
        {momFactors}];
    Catch[
      {
        "slash",
        Sort[
          Function[trm,
              Module[{q, c},
                q = Cases[trm, ntVec[qq_, m_] /; m === mu :> qq, {0, Infinity}];
                c = trm /. ntVec[_, m_] /; m === mu :> 1;
                If[Length[q] =!= 1 || !FreeQ[c, ntVec] || !NumericQ[c],
                  Throw[$Failed, "sigleg"]];
                {c, First[q]}]
            ] /@ terms]},
      "sigleg"]];

(* chain-order the 2 gammas of one commutator term (spinor in->out) → {legFirst, legSecond, din, dout,
   scalar}, or $Failed if the term is not a clean isolated 2-γ pair. *)

ntSigmaTermInfo[term_] := Module[{tf, gs, g1, g2, first, second, din, dout},
    tf =
      If[Head[term] === Times,
        List @@ term,
        {term}];
    gs = Cases[tf, _ntGamma];
    If[Length[gs] =!= 2,
      Return[$Failed]];
(* the two γ's (and their ntVecs) must be the term's ONLY tensor structure — no colour, projector,
   metric, other Dirac heads: anything else means the Plus is not a clean bare commutator. *)
    If[!FreeQ[tf, _ntGamma5 | _ntDeltaDirac | _ntSigma | _ntSUNT | _ntSUNDeltaFund | _ntSUNf | _ntSUNDeltaAdj | _ntTransProj | _ntLongProj | _ntMetric | _ntEpsilon],
      Return[$Failed]];
    {g1, g2} = gs;
    Which[
      g1[[3]] === g2[[2]],
        first = g1;
        second = g2;
        din = g1[[2]];
        dout = g2[[3]],
      g2[[3]] === g1[[2]],
        first = g2;
        second = g1;
        din = g2[[2]];
        dout = g1[[3]],
      True,
        Return[$Failed]];
    {ntSigmaLeg[first, tf], ntSigmaLeg[second, tf], din, dout, Times @@ Select[tf, scalarQ]}];

(* {scalar, ntSigma[...]} for a 2-term Plus that is a bare γ-commutator, else $Failed. The two terms
   must share spinor endpoints, carry the SAME two legs in SWAPPED order, and have opposite scalar.
   Result uses term-1's leg order and scalar: c·([X,Y] block) (sign-symmetric — either term gives the
   same fold, since c_2 = −c_1 and swapping legs negates the bracket). *)

ntRecognizeComm[p_] := Module[{terms, a, b},
    If[Head[p] =!= Plus || Length[p] =!= 2,
      Return[$Failed]];
    terms = List @@ p;
    a = ntSigmaTermInfo[terms[[1]]];
    b = ntSigmaTermInfo[terms[[2]]];
    If[a === $Failed || b === $Failed,
      Return[$Failed]];
    If[a[[3]] =!= b[[3]] || a[[4]] =!= b[[4]],
      Return[$Failed]
    ];(* same spinor endpoints *)
    If[!(a[[1]] === b[[2]] && a[[2]] === b[[1]]),
      Return[$Failed]
    ];(* legs swapped between terms *)
    If[a[[1]] === a[[2]],
      Return[$Failed]
    ];(* identical legs ⇒ [X,X]=0 *)
    If[Simplify[a[[5]] + b[[5]]] =!= 0,
      Return[$Failed]
    ];(* opposite scalar sign *)
    {a[[5]], ntSigma[a[[1]], a[[2]], a[[3]], a[[4]]]}];

(* replace every recognized bare γ-commutator Plus factor with its single ntSigma token + scalar.
   Set env NT_NO_SIGMA_FOLD=1 to disable the optimization (the commutator then distributes into two
   traces as before) — a safety switch and the baseline for measuring the fold's sub-term reduction. *)

$ntSigmaFold = (Environment["NT_NO_SIGMA_FOLD"] === $Failed);

(* Fold every bare γ-commutator Plus among `factors` into a single ntSigma token (gated by
   $ntSigmaFold), recursively: find a factor that is a Plus recognisable as a commutator [A,B]
   (ntRecognizeComm), replace it with its two σ legs, recurse until none remain — so the antisymmetric
   γ-pair is never distributed into two separate Dirac traces. *)

foldDiracSigma[factors_List] := Module[{commPlus, recognized},
    If[!$ntSigmaFold,
      Return[factors]];
    commPlus = FirstCase[factors, q_ /; (Head[q] === Plus && ntRecognizeComm[q] =!= $Failed), Missing[], {1}];
    If[MissingQ[commPlus],
      Return[factors]];
    recognized = ntRecognizeComm[commPlus];
    foldDiracSigma[Join[DeleteCases[factors, commPlus, {1}, 1], {recognized[[1]], recognized[[2]]}]]];

(* Split a colour-ENTANGLED Lorentz/Dirac component (the AqbqDirect147 4/7 quark-gluon vertex, where
   colour and Dirac do NOT factor globally, only per term) into per-colour-structure groups so the
   Lorentz-only builderInv can lower each. Algorithm:
     1. fold σ commutators, then EXPAND only the entangled Pluses (those mixing colour with
        Dirac/Lorentz) into branches — single-sector sums stay eager (no monomial blow-up).
     2. each branch → {colourProduct, {core, scal, restStr}} via compileDirac / chunkLorInv.
     3. GatherBy colour product, emitting ONE {colourNet, bodyNets, scalar, restNets} entry per group
        (so the net count tracks the colour graph, not the branch×ordering explosion); the generator
        sums each group's branches at runtime. *)

splitColourGroupsInv[factors0_, ids_, env_, mask_, nc_] :=
  Module[{factors = foldDiracSigma[factors0], entangledQ, needExpand, keepAll, distributed, terms, branchNets, groups},
    entangledQ[x_] := Head[x] === Plus && (colourEntangledQ[x] || !FreeQ[x, _ntGamma | _ntGamma5 | _ntDeltaDirac]);
    needExpand = Select[factors, entangledQ];
    keepAll = Select[factors, !entangledQ[#]&];
    distributed = Expand[Times @@ needExpand];(* small: product of the entangled Pluses only *)
    terms =
      If[Head[distributed] === Plus,
        List @@ distributed,
        {distributed}];
    (* per branch -> {colourProduct, list-of-{bodyOrCore,scal,restStr}} *)
    branchNets =
      Function[term,
          Module[{termFactors, colProd, rest},
            termFactors =
              Join[
                If[Head[term] === Times,
                  List @@ term,
                  {term}],
                keepAll];
            colProd = Times @@ Cases[termFactors, ctFacInv];
            rest = DeleteCases[termFactors, ctFacInv];(* gammas + Lorentz + numeric coeff (no colour) *)
(* Level-1 DeleteCases only strips BARE group-head factors. Anything that buries one (a Power, or
   a Plus that entangledQ did not expand) leaves colour in the remainder, which then reaches the
   Lorentz-only builderInv and is CForm'd into the .cpp. Fail here, where the offender is still
   identifiable, rather than at the emission chokepoint three layers away. *)
            If[!FreeQ[rest, Alternatives @@ ctHeadsInv],
              Message[MakeNTKernel::colrest, Short[DeleteDuplicates @ Cases[rest, Alternatives @@ ctHeadsInv, {0, Infinity}], 6], Short[rest, 8]];
              Abort[]];
            {
              colProd,
              If[!FreeQ[rest, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac | _ntDressedNum | _ntDiracSlot],
                {compileDirac[rest, ids, env, mask, nc]},
                (* gamma chain: {core, scal, projectorRest} *)
                ({#[[1]], #[[2]], ""}&) /@ chunkLorInv[Times @@ rest, ids, env, mask, nc]]}]
        ] /@ terms;
    groups = GatherBy[branchNets, First];
(* one {colourNet, bodyNet, scalar, restNet} entry per colour group.
   INV backend: fuse a group's branches into ONE NetVal body via add/sc/contract — valid because
   each branch's `core` is itself a NetVal (dirac_value(...)); restNet is a string.
   NUMERIC backend: a branch's `core` is a BARE `DiracNet{...}` (contracted against its Lorentz
   rest only at runtime by numeric_value_netval), so fusing it into `add(sc(DiracNet{...}), …)` is
   a C++ type error (the `1/4/7` bug). Keep ONE entry per colour group (so the net count tracks the
   graph, not the branch×ordering explosion), but carry the group's branches as a LIST of sub-terms:
   bodyNet = {core_b…} (each a DiracNet literal, or a Lorentz NetVal for a gamma-free branch),
   restNet = {{rest_b, scal_b}…} (parallel). The generator sums mp[i]=Σ_b scal_b·numeric_value_netval(
   dnet_b, lnet_b) — matching how the inv body sums its branches, but type-correct for numeric. *)
    Function[group,
        Module[{colProd = group[[1, 1]], colNet, colScalar, branchRecs},
          {colNet, colScalar} =
            If[colProd === 1,
              {"SUNNet{}", 1},
              compileColG[colProd, ids]];
          branchRecs = Flatten[group[[All, 2]], 1];(* each = {core, scal, restStr} *)
          {colNet, branchRecs[[All, 1]], colScalar, ({#[[3]], #[[2]]}&) /@ branchRecs}]
      ] /@ groups];

(* ---- colour/group dialect: a constant SU(N) component is a product of structure constants,
        generators, and Kronecker deltas; emit it as a `SUNNet` literal for the generator's
        numeric SU(N) contraction (the et engine can't instantiate the four-gluon colour tensor
        type). Each head carries its own group rank N as the leading argument, so one net can mix
        several groups (e.g. colour SU(Nc) ⊗ flavour SU(Nf)) — sun_net.hpp's sun_value_cx
        contracts each rank separately and multiplies. Returns {colourNetString, factoredScalar}. *)
(* Each factor is minted by the per-rank `sun<n>` SUNEnv (declared in the generator main, one per
   distinct rank appearing in the colour nets), so the group rank is written once, not on every factor.
   `sun<n>` is the network::SUNEnv analogue of the numeric LorentzEnv. *)

colFacG[ntSUNf[n_, a_, b_, c_], ids_] :=
  "sun" <> ToString[n] <> ".f(" <> ToString[ids[a]] <> "," <> ToString[ids[b]] <> "," <> ToString[ids[c]] <> ")";

colFacG[ntSUNDeltaAdj[n_, a_, b_], ids_] :=
  "sun" <> ToString[n] <> ".deltaAdj(" <> ToString[ids[a]] <> "," <> ToString[ids[b]] <> ")";

colFacG[ntSUNT[n_, a_, i_, j_], ids_] :=
  "sun" <> ToString[n] <> ".T(" <> ToString[ids[a]] <> "," <> ToString[ids[i]] <> "," <> ToString[ids[j]] <> ")";

colFacG[ntSUNDeltaFund[n_, i_, j_], ids_] :=
  "sun" <> ToString[n] <> ".deltaFund(" <> ToString[ids[i]] <> "," <> ToString[ids[j]] <> ")";

(* per-component diagonal dressings: parse the spec into a per-component dressing-id vector
   (component → dr, -1 = drop; 1-based physics indices) and emit a diag factor carrying it. *)

colFacG[ntSUNDiagFund[n_, i_, j_, spec_, scale_], ids_] :=
  "sun" <> ToString[n] <> ".diagFund(" <> ToString[ids[i]] <> "," <> ToString[ids[j]] <> "," <> diagVecStr[diagComp2Dr[spec, scale, n]] <> ")";

colFacG[ntSUNDiagAdj[n_, a_, b_, spec_, scale_], ids_] :=
  "sun" <> ToString[n] <> ".diagAdj(" <> ToString[ids[a]] <> "," <> ToString[ids[b]] <> "," <> diagVecStr[diagComp2Dr[spec, scale, n^2 - 1]] <> ")";

(* CATCH-ALL, and it must stay LAST: the six rules above are the only lowerable colour factors.
   Without this a non-matching factor returns UNEVALUATED and StringRiffle happily ToString's it
   into the generator .cpp — `colFacG[Plus[...], <|v1 -> 0, ...|>]` — which is how the four-quark
   Fierz flavour sum surfaced (as a clang syntax error, three layers away from its cause).
   compileTInv has had exactly this guard (MakeNTKernel::tleak) since the ZAAqbq metric leak;
   colFacG was the one per-head dispatcher in this file with neither a Plus branch nor a
   fallthrough. Fail loudly and locally instead. *)

colFacG[e_, _] := (
    Message[MakeNTKernel::colleak, Short[e, 6]];
    Abort[]);

compileColG[e_, ids_] := Module[
    {
      parts =
        If[Head[e] === Times,
          List @@ e,
          {e}],
      sc,
      tn},
(* A colour/flavour SUM raised to a power cannot be expanded by repetition: the copies would share
   index labels (see MakeNTKernel::colpow). Refuse before the rewrite below can do it. *)
    Cases[
      parts,
      Power[b_Plus, k_Integer?Positive] /; !scalarQ[b] :>
        (
          Message[MakeNTKernel::colpow, k, Short[b, 6]];
          Abort[])];
(* a colour/flavour factor raised to an integer power (e.g. deltaAdjFlav^2 from a CLOSED meson
   loop: delta_adj(a,b)^2 -> the flavour trace N^2-1) must be expanded into repeated SUNNet
   factors so the C++ sun_value contracts the shared indices — colFacG handles a single head,
   not Power[head,k], so an un-expanded power leaks Mathematica syntax into the generator.
   Restricted to a BARE group head: the old `! scalarQ[b]` guard also admitted Power[Plus[..],k]
   and Power[Times[..],k], whose repetition duplicates labels rather than closing a self-trace. *)
    parts = parts /. Power[b_, k_Integer?Positive] /; MatchQ[b, Alternatives @@ ctHeadsInv] :> Sequence @@ ConstantArray[b, k];
    sc = Select[parts, scalarQ];
    tn = Select[parts, !scalarQ[#]&];
    {"SUNNet{" <> StringRiffle[colFacG[#, ids]& /@ tn, ", "] <> "}", Times @@ sc}];

(* ---- a colour/flavour component that is a SUM ------------------------------------------------
   `SUNNet` is a flat PRODUCT (std::vector<SUNFac>) with no sum node and no coefficient field, and
   mergeColNet splices `SUNNet{...}` literals by string surgery — so a sum cannot be represented
   inside one net, and there is no colour analogue of compileTInv's `add(...)`.

   It does not need one. Colour folds to a SCALAR (sun_value_cx -> Cx), and the generator already
   sums colour structures by emitting SEVERAL nets that share a group:
       for(int d: grp) acc = acc + mp[d]*env.constant(colv[d]);
   So a summed colour component lowers to a LIST of {net, scalar} branches — the same shape
   chunkLorInv already returns and its callers already loop over.

   ONE Expand over the product of all the diagram's constant components does the cross-product in
   one step (the same device splitColourGroupsInv uses for its entangled sums), so several summed
   constant components do not need pairwise outer products.

   Reached only from the CONSTANT-component path: splitColourGroupsInv's colProd is a level-1
   Cases over an already-Expanded branch, hence flat by construction, and keeps using compileColG. *)

$ntColSumMaxBranches = 4096;

MakeNTKernel::colsum = "compileColGSum: a constant colour/flavour component expands to `1` summed branches (limit `2`). Each branch costs one emitted SUNNet and one net record, so this would blow up the generator. Raise $ntColSumMaxBranches if the flow genuinely needs it.";

compileColGSum[e_, ids_] := Module[{
    terms =
      With[{x = Expand[e]},
        If[Head[x] === Plus,
          List @@ x,
          {x}]]},
    If[Length[terms] > $ntColSumMaxBranches,
      Message[MakeNTKernel::colsum, Length[terms], $ntColSumMaxBranches];
      Abort[]];
    compileColG[#, ids]& /@ terms];

(* ---- et-vs-numeric size guard for a FUNDAMENTAL colour/flavour component ----------
   A small generator contraction (a quark loop's T^a) instantiates fine as an et ETensor and folds
   to a constexpr; a large one OOMs the compiler exactly like the four-gluon adjoint type. The metric
   below is the widest intermediate of contract_all under the greedy order — the product of the
   still-open axis extents (adjoint N^2-1, fundamental N) at each step — which directly predicts the
   ETensor type size. Above the threshold the component takes the numeric SUNNet path instead. *)

$NumTracerFundETMaxEntries = 200000;

(* `tensorQ` is HEAD-matching with a `tensorQ[_] = False` catch-all, so tensorQ[Plus[...]] is False —
   a Plus is neither scalarQ nor tensorQ. `Select[facs, tensorQ]` therefore SILENTLY DROPS a sum
   vertex, which here would understate the contraction width and could route a component to the et
   path that belongs on the numeric one (the OOM $NumTracerFundETMaxEntries exists to prevent).
   Flatten a Plus to its widest summand instead: every summand shares the sum's free indices, so the
   summand with the most/widest labels bounds the whole sum. *)

tensorPartsOf[facs_] := Flatten[
    Function[f,
        Which[
          tensorQ[f],
            {f},
          Head[f] === Plus,
            tensorPartsOf[List @@ MaximalBy[List @@ f, LeafCount][[1]]],
          Head[f] === Times,
            tensorPartsOf[List @@ f],
          Head[f] === Power && IntegerQ[f[[2]]] && f[[2]] >= 1 && !scalarQ[f],
            ConstantArray[tensorPartsOf[{f[[1]]}], f[[2]]],
          True,
            {}]
      ] /@ facs];

labelDimAssoc[facs_] := Module[{assoc = <||>},
    Function[h,
        With[{n = sunRankOf[h]},
          Switch[h,
            _ntSUNf | _ntSUNDeltaAdj | _ntSUNDiagAdj,
              (assoc[#] = n^2 - 1)& /@ labelsOf[h],
            _ntSUNT,
              (
                assoc[First[labelsOf[h]]] = n^2 - 1;
                (assoc[#] = n)& /@ Rest[labelsOf[h]]),
            _ntSUNDeltaFund | _ntSUNDiagFund,
              (assoc[#] = n)& /@ labelsOf[h]]]
      ] /@ tensorPartsOf[facs];
    assoc];

fundMetric[facs_] :=
  Module[{tn = orderFactors[tensorPartsOf[facs]], dimA = labelDimAssoc[facs], cnt = <||>, widest = 1, open},
    Do[
      (cnt[#] = Lookup[cnt, #, 0] + 1)& /@ labelsOf[f];
      open = Keys[Select[cnt, OddQ]];
      widest = Max[widest, Times @@ (dimA /@ open)],
      {f, tn}];
    widest];

(* ---- Dirac trace -> Lorentz NetVal (for the invariant-basis path) ----------
   The et engine takes a spinor trace by closing a cyclic chain; the inv generator works in
   scalar products, so we instead apply the gamma-trace THEOREM symbolically, turning the
   closed chain of gammas into a signed sum of metric products. Each slashed momentum is
   ntGamma[mu,..] * ntVec[q,mu], so the metric mu-pairings contract the vecs/projectors that
   share those Lorentz labels — leaving a pure-Lorentz network the existing reduce->rebase
   pipeline handles unchanged. (No gamma5 support yet; gamma5-bearing traces are not expanded.) *)

diracIn[ntGamma[_, a_, _]] := a; diracOut[ntGamma[_, _, b_]] := b;

diracIn[ntGamma5[a_, _]] := a; diracOut[ntGamma5[_, b_]] := b;

diracIn[ntSigma[_, _, a_, _]] := a; diracOut[ntSigma[_, _, _, b_]] := b;

diracIn[ntDeltaDirac[a_, _]] := a; diracOut[ntDeltaDirac[_, b_]] := b;

diracIn[ntDressedNum[_, a_, _]] := a; diracOut[ntDressedNum[_, _, b_]] := b;

diracIn[ntDiracSlot[_, a_, _, _]] := a; diracOut[ntDiracSlot[_, _, b_, _]] := b;

(* walk the closed spinor loop, returning trace-ordered tokens: a γ's Lorentz label, or the marker
   $g5 for a γ5 (a spinor-δ just connects indices and contributes no token). *)

orderDiracChain[facs_] := Module[{byIn, cur, toks = {}, seen = {}},
    byIn = Association[(diracIn[#] -> #)& /@ facs];
    cur = First[facs];
    While[
      !MemberQ[seen, cur],
      AppendTo[seen, cur];
      Which[
        MatchQ[cur, _ntGamma],
          AppendTo[toks, First[cur]],
        MatchQ[cur, _ntGamma5],
          AppendTo[toks, $g5]];
      cur = Lookup[byIn, diracOut[cur], Missing[]];
      If[MissingQ[cur],
        Break[]]];
    toks];

(* tr[g^m1 ... g^m2n] / 4 = sum over pairings of signed products of metrics (recursive). The pairing
   STRUCTURE depends only on the chain LENGTH, not the labels: gammaTraceSum[{a,b,c,d}] is
   gammaTraceSum[{1,2,3,4}] with positions -> labels. A quark triangle distributes into hundreds of
   branches that each trace a long (10-12 gamma) chain — recomputing the (2n-1)!! pairing sum per
   branch was catastrophic (minutes / GBs per component). Instead compute the canonical positional
   trace ONCE per length (memoized) and substitute each branch's labels — hundreds of cheap
   substitutions replace hundreds of exponential recomputations. *)

rawGammaTrace[{}] := 1;

rawGammaTrace[lst_] :=
  Sum[(-1) ^ j ntMetric[lst[[1]], lst[[j]]] rawGammaTrace[Delete[lst, {{1}, {j}}]], {j, 2, Length[lst]}];

canonGammaTrace[n_] := canonGammaTrace[n] =
  rawGammaTrace[Range[n]];(* positions 1..n, memoized per length *)

gammaTraceSum[{}] := 1;

gammaTraceSum[lst_] := canonGammaTrace[Length[lst]] /. Dispatch[Thread[Range[Length[lst]] -> lst]];

(* tr[γ5 g^m1 ... g^m2n] / 4 = the chiral trace. Vanishes for < 4 gammas; the 4-gamma base is the
   Levi-Civita ε — the normalisation tr(γ5 γ γ γ γ) = 4 ε is pinned against the typed-out component
   engine (test_inv_eps, constant C = 4). ≥ 6 gammas peel one metric exactly like the ordinary trace,
   with the ε base seeding the recursion (the standard FORM chiral reduction). *)

gammaTraceSum5[lst_] := Which[
    Length[lst] < 4,
      0,
    Length[lst] == 4,
      ntEpsilon[lst[[1]], lst[[2]], lst[[3]], lst[[4]]],
    True,
      Sum[(-1) ^ j ntMetric[lst[[1]], lst[[j]]] gammaTraceSum5[Delete[lst, {{1}, {j}}]], {j, 2, Length[lst]}]];

(* trace of a token chain: collapse any number of γ5 to ≤ 1 (γ5 γ^μ = −γ^μ γ5, γ5² = 1). The sign to
   slide every γ5 together past the intervening γ's is (−1)^(#{γ before each γ5}); an even γ5 count →
   the ordinary trace, an odd count → the chiral trace. (This is what guarantees the C++ ε-network's
   "≤ 1 ε per term" invariant.) *)

diracTrace[toks_] := Module[{n5 = Count[toks, $g5], gammas = DeleteCases[toks, $g5], sign},
    sign = (-1) ^ Total[(Length[Cases[Take[toks, # - 1], Except[$g5]]])& /@ Flatten[Position[toks, $g5]]];
    sign *
      If[EvenQ[n5],
        gammaTraceSum[gammas],
        gammaTraceSum5[gammas]]];

(* expand a Dirac-trace component (factor list) to a pure-Lorentz expression; non-Dirac
   components pass through unchanged. The tr(1)=4 is the explicit leading factor.
   LEGACY: the numeric backend traces in C++ (compileDirac → network::dirac_value); this symbolic
   path is retained only for cross-validation. *)

expandDiracComponent[factors_List] := Module[{dir, rest},
    dir = Select[factors, MatchQ[#, _ntGamma | _ntGamma5 | _ntDeltaDirac]&];
    rest = DeleteCases[factors, _ntGamma | _ntGamma5 | _ntDeltaDirac];
    If[dir === {},
      Times @@ factors,
      4 diracTrace[orderDiracChain[dir]] (Times @@ rest)]];

(* ---- Dirac trace in the C++ generator (et/inv/dirac.hpp `dirac_value`) -----------------------
   The inv backend emits the gamma chain as a `DiracNet` literal and the generator contracts the
   closed spinor loop NUMERICALLY into a Lorentz net (metrics/vectors over the free gluon legs), then
   `reduce` folds it with the projector — instead of expanding the (2n−1)!! pairing sum symbolically
   in Mathematica (which exploded for the 3-gluon-vertex quark triangle). This parallels how colour is
   folded by `SUNNet`/`compileColG`. *)

$ntDiracFree = 900000;(* fresh Lorentz-label base for slash–slash pairings: above every component id *)
(* the gamma/gamma5 factors in trace order (walk the closed spinor loop; spinor-δ connectors carry no
   token but are followed). Like `orderDiracChain` but returns the factors so each can be classified.

   UNDIRECTED cycle walk: each Dirac factor is an EDGE between its two spinor labels (spinorLabelsHead);
   every label in a closed spinor loop has degree 2, so the walk is deterministic. We seed at
   `First[facs]` ENTERING on its `diracIn` label (so it leaves on `diracOut`) and at each step leave the
   current factor by its OTHER endpoint, picking the unique unvisited neighbour there. For a
   consistently-oriented loop (every node one in / one out — Zq/ZA/ZA3/ZA4/ZAqbq, which close via an
   oriented γ) this reproduces the old directed `byIn` forward walk token-for-token (byte-identical).
   The robustness is for a loop closed by a SYMMETRIC external spinor-δ (`ntDeltaDirac[d1,d2]`, e.g. the
   σL scalar external projector): the old directed walk collided on the shared `diracIn`/`diracOut` keys
   and could not traverse the cycle, collapsing the trace to a single γ5 → tr(γ5)=0. The undirected
   adjacency closes the cycle correctly regardless of start/orientation (the trace is cyclic). *)

orderDiracFacs[facs_] :=
  Module[{nodeFacs = Association[], cur = 1, prevLabel, out = {}, seen = {}, labels, exitLabel, nexts, nTok},
    Do[
      Module[{ls = spinorLabelsHead[facs[[i]]]},
        (nodeFacs[#] = Append[Lookup[nodeFacs, #, {}], i])& /@ ls],
      {i, Length[facs]}];
    prevLabel = diracIn[facs[[1]]];(* enter First on its diracIn so it exits on diracOut *)
    While[
      !MemberQ[seen, cur],
      AppendTo[seen, cur];
      If[MatchQ[facs[[cur]], _ntGamma | _ntGamma5 | _ntSigma | _ntDressedNum | _ntDiracSlot],
        AppendTo[out, facs[[cur]]]];
      labels = spinorLabelsHead[facs[[cur]]];
      exitLabel = First[DeleteCases[labels, prevLabel], Missing[]];(* the OTHER endpoint *)
      If[MissingQ[exitLabel],
        Break[]];
      nexts = Select[DeleteCases[Lookup[nodeFacs, exitLabel, {}], cur], !MemberQ[seen, #]&];
      If[nexts === {},
        Break[]];
      prevLabel = exitLabel;
      cur = First[nexts]];
(* LOUD GUARD: a closed spinor loop must consume EVERY token-bearing factor. If the walk fell short
   (a fragmented / improperly-closed loop), the emitted trace would silently drop γ structure (the
   historical hSigL meson-sector collapse: a single surviving γ5 → tr(γ5)=0). Abort rather than emit a
   wrong kernel. The δ-only "connector" factors carry no token, so they are excluded from the count. *)
    nTok = Count[facs, _ntGamma | _ntGamma5 | _ntSigma | _ntDressedNum | _ntDiracSlot];
    If[Length[out] =!= nTok,
      Print["[NumTracer] ERROR: orderDiracFacs walk consumed ", Length[out], " of ", nTok, " token-bearing Dirac factors — a spinor loop did not close (would silently drop γ ", "structure → collapsed trace). Loop factors:\n  ", facs];
      Abort[]];
    out];

(* A component may contain SEVERAL independent closed spinor loops (a quark loop + the
   projection-closed external line, tied together only by gluon propagators). orderDiracFacs walks ONE
   loop and would drop the rest; partition the Dirac factors into spinor-connected groups first
   (two factors share a loop iff they share a spinor index), then order each loop. Returns a LIST of
   ordered token lists (one per loop) — a single-loop component yields a one-element list. The C++
   contraction traces each loop separately and contracts their shared gluon legs (DFac::LoopSep). *)
(* MEMOISED: called once per Dirac component, but the distinct Dirac structures are far fewer than the
   calls — the same chains recur across diagrams and colour branches (the same redundancy the net-term
   CSE exploits downstream). Building the spinor-loop graph on every call, rather than once per distinct
   chain, dominated the net-build on a dense dressed flow. Pure function of `facs`; $odCache is cleared
   per generation in mkGenerateKernel. *)

$odCache = <||>;

orderDiracLoops[facs_] := With[{h = Hash[facs]},
    Lookup[$odCache, h, $odCache[h] = orderDiracLoopsBody[facs]]];

orderDiracLoopsBody[facs_] := If[Length[facs] <= 1,
    {orderDiracFacs[facs]},
    Module[
      {sp = spinorLabelsHead /@ facs, edges, g, comps},
(* Two factors are adjacent iff they SHARE a spinor label, so index label -> nodes and read the edges
   off the buckets: O(n) in the chain length. Testing all Subsets[...,{2}] pairs instead is O(n^2),
   which was the steepest term in the net-build and the one that would bite hardest on flows with
   longer Dirac chains. Sort+DeleteDuplicates reproduce the old pair scan's output exactly —
   lexicographic order (so the Graph, hence ConnectedComponents' component order, is unchanged) and
   each pair once (two factors may share BOTH labels). *)
      edges = Sort @ DeleteDuplicates @ Flatten[Subsets[#, {2}]& /@ Values @ GroupBy[Flatten[Table[{l, i}, {i, Length[facs]}, {l, sp[[i]]}], 1], First -> Last], 1];
      g = Graph[Range[Length[facs]], UndirectedEdge @@@ edges];
      comps = ConnectedComponents[g];
(* a SINGLE spinor loop → the exact original path (orderDiracFacs[facs]), so existing flows stay
   byte-identical; only genuinely multi-loop components are split (indices sorted to preserve the
   original relative order within each loop). *)
      If[Length[comps] <= 1,
        {orderDiracFacs[facs]},
        (orderDiracFacs[facs[[Sort[#]]]])& /@ comps]]];

(* {netString, scalar} for a colour-free component (gammas + Lorentz + numeric coeff). The gamma chain
   becomes `dirac_value(DiracNet{...})` (each γ: a SLASH `dslash` if its Lorentz index is contracted
   with an `ntVec[q,μ]`, else a FREE leg `dgamma(μ)` that contracts the projector); the remaining
   (Lorentz) factors compile through `compileTInv` and are `contract`ed with the trace. A component
   with no Dirac factor is just its Lorentz net. *)
(* frame resolver for dressed-numerator option coefficients (ntSP/ntVec[q,i] -> components). Set in
   mkGenerateKernel to the diagram's resolveScale; Identity when the dressed path is inactive. *)

$ntDressResolve = Identity;

(* the C++ DSlot literal for one ntDressedNum: a sum of DSlotOpt{Cx coeff, {dress ids}, slash?, vlc}.
   Each option's scalar coefficient is frame-resolved then split into a complex number × dressing atoms
   (drDecompose); a "slash" structure's momenta become a vlc of {1.0, env Base} pairs (coeff-1, like the
   single-vec slash above), an "ident" structure is the spinor identity (slash=false, empty vlc). *)
(* MEMOISED: called once per dressed-numerator token, so a dense dressed flow makes tens of thousands of
   calls — but has only a handful of DISTINCT dressed numerators, since the same propagator numerator
   recurs in every diagram and colour branch. Resolving/decomposing one is expensive, and unmemoised this
   was the single largest cost in the net-build. Pure given `env` and $ntDressResolve, both fixed for the
   generation; $dsCache is cleared per generation in mkGenerateKernel alongside $ctCache. *)

$dsCache = <||>;

dressedSlotStr[gf : ntDressedNum[_, _, _], env_] := With[{h = Hash[{gf, env}]},
    Lookup[$dsCache, h, $dsCache[h] = dressedSlotStrBody[gf, env]]];

(* Returns the LIST of per-option "DSlotOpt{…}" strings (NOT the wrapped "DSlot{…}"). The generator
   (emitNumericGenerator) expands the Cartesian product of the chain's slots' options into one
   single-option sub-term per combination, so each slot is delivered here as its bare option list. *)
dressedSlotStrBody[ntDressedNum[opts_, _, _], env_] := Function[opt,
          Module[{num, dr, vlcStr},
            {num, dr} = drDecompose[$ntDressResolve[opt[[1]]]];
            vlcStr =
              If[opt[[2, 1]] === "slash",
                StringRiffle[
                  Function[cq,
                      Module[{c = cq[[1]], q = cq[[2]]},
                        If[!KeyExistsQ[env, q],
                          Print["[NumTracer] ERROR: dressed slash momentum ", q, " absent from env ", Keys[env]];
                          Abort[]];
                        "{" <> cppNum[c] <> "," <> ToString[env[q]["Base"]] <> "}"]
                    ] /@ opt[[2, 2]],
                  ", "],
                ""];
            (* LEVER (b): return the STRUCTURE separately from the dressing, as {structStr, num, dr}.
                     structStr is the dressing-free DSlotOpt (coeff 1, no dress atoms) — ident → empty toks,
                     slash → one dslash token (netFacs empty: a k=0 propagator numerator has no open leg).
                     num (numeric Cx) folds into the sub-term scalar and dr (dress-atom ids) becomes the
                     DPoly key, so the trace table dedups on structure alone. *)
            {"DSlotOpt{Cx{1,0}, {}, {" <>
              If[opt[[2, 1]] === "slash", "dslash({" <> vlcStr <> "})", ""] <> "}, {}}", num, dr}]
        ] /@ opts;

(* ---- general collected Dirac slot → C++ DSlot literal (Stage 4, any open-leg count) ------------
   An ntDiracSlot option keeps its WHOLE structure (Dirac chain × Lorentz-net factors); here we split
   each option into DSlotOpt{coeff, {dress}, {toks}, {netFacs}} — the Dirac chain as a token list
   (dgamma/dslash/dcomm/dg5, open legs = ids of the free Lorentz tokens) and the Lorentz factors as
   network::Elem literals (the gluon propagator/metric that closes the open leg). Every label already
   has an id and every momentum an env slot (allLabels/momentumOf recurse into the slot via
   Cases[Infinity]); internal legs (the γ↔projector bridge) keep their own distinct ids, closed within
   the option, so no fresh-id allocation is needed. *)
(* Order one OPEN spinor chain din→dout, returning the token-bearing factors (γ/γ5/σ) in chain order;
   δ connectors are followed but carry no token. Walk from the din endpoint (degree 1) along spinor
   adjacency — like orderDiracFacs but seeded at a KNOWN endpoint (the chain is open, not a cycle). *)
orderOpenChain[facs_, din_] :=
  Module[{nodeFacs = <||>, cur, prevLabel = din, out = {}, seen = {}, exitLabel, nexts, start},
    Do[(nodeFacs[#] = Append[Lookup[nodeFacs, #, {}], i]) & /@ spinorLabelsHead[facs[[i]]], {i, Length[facs]}];
    start = Lookup[nodeFacs, din, {}];
    If[start === {}, Return[Select[facs, MatchQ[#, _ntGamma | _ntGamma5 | _ntSigma] &]]]; (* fallback: no din endpoint *)
    cur = First[start];
    While[! MemberQ[seen, cur],
      AppendTo[seen, cur];
      If[MatchQ[facs[[cur]], _ntGamma | _ntGamma5 | _ntSigma], AppendTo[out, facs[[cur]]]];
      exitLabel = First[DeleteCases[spinorLabelsHead[facs[[cur]]], prevLabel], Missing[]];
      If[MissingQ[exitLabel], Break[]];
      nexts = Select[DeleteCases[Lookup[nodeFacs, exitLabel, {}], cur], ! MemberQ[seen, #] &];
      If[nexts === {}, Break[]];
      prevLabel = exitLabel; cur = First[nexts]];
    out];

$dslCache = <||>;
diracSlotStr[gf : ntDiracSlot[_, _, _, _], ids_, env_, mask_, nc_] := With[{h = Hash[{gf, ids, env}]},
    Lookup[$dslCache, h, $dslCache[h] = diracSlotStrBody[gf, ids, env, mask, nc]]];

(* Returns the LIST of per-option "DSlotOpt{…}" strings (see dressedSlotStrBody): the generator
   expands the Cartesian product of a chain's slots into single-option sub-terms. *)
diracSlotStrBody[ntDiracSlot[opts_, din_, dout_, legs_], ids_, env_, mask_, nc_] := Function[opt,
        Module[{num, dr, facs, vecOf, gammaLegs, diracFacs, lorFacs, ordered, toks, netFacs, legStr, sigStr},
          {num, dr} = drDecompose[$ntDressResolve[opt[[1]]]];
          facs = If[Head[opt[[2]]] === Times, List @@ opt[[2]], {opt[[2]]}];
          (* μ -> slash momentum q (an ntVec sharing a γ's Lorentz leg) *)
          vecOf = Association[Reverse[Cases[facs, ntVec[q_, m_] :> (m -> q)]]];
          gammaLegs = Cases[facs, ntGamma[gm_, _, _] :> gm];
          diracFacs = Select[facs, MatchQ[#, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac] &];
          (* Lorentz-net factors = non-Dirac tensors that are NOT a slash-vec (a slash-vec's μ is a γ leg) *)
          lorFacs = Select[facs, (tensorQ[#] && ! MatchQ[#, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac] &&
                          ! MatchQ[#, ntVec[_, m_ /; MemberQ[gammaLegs, m]]]) &];
          (* an ntSigma leg → C++ arg (free open leg id, or a slashed-leg vlc) *)
          legStr[{"slash", pairs_List}] :=
            "{" <> StringRiffle[("{" <> cppNum[#[[1]]] <> "," <> ToString[env[#[[2]]]["Base"]] <> "}") & /@ pairs, ", "] <> "}";
          legStr[{"free", mu_}] := ToString[ids[mu]];
          sigStr[legA_, legB_] := Module[{ta = legA[[1]], tb = legB[[1]], sa = legStr[legA], sb = legStr[legB]},
            Which[ta === "free" && tb === "free", "dcomm(" <> sa <> ", " <> sb <> ")",
                ta === "slash" && tb === "slash", "dcomm_ss(" <> sa <> ", " <> sb <> ")",
                ta === "free" && tb === "slash", "dcomm_fs(" <> sa <> ", " <> sb <> ")",
                True, "dcomm_sf(" <> sa <> ", " <> sb <> ")"]];
          ordered = orderOpenChain[diracFacs, din];
          toks = Function[gf2,
            Which[
              MatchQ[gf2, _ntGamma5], "dg5()",
              MatchQ[gf2, _ntSigma], sigStr[gf2[[1]], gf2[[2]]],
              True, With[{mu = First[gf2]},   (* ntGamma: slash if its leg is an ntVec momentum, else a free open leg *)
                If[KeyExistsQ[vecOf, mu],
                  "dslash({{1.0," <> ToString[env[vecOf[mu]]["Base"]] <> "}})",
                  "dgamma(" <> ToString[ids[mu]] <> ")"]]]] /@ ordered;
          netFacs = builderInvElem[#, ids, env] & /@ lorFacs;
          (* LEVER (b): {structStr, num, dr} — dressing-free DSlotOpt (coeff 1, no dress) + the numeric
                  Cx (num) and dress-atom ids (dr) carried separately by the sub-term. See dressedSlotStrBody. *)
          {"DSlotOpt{Cx{1,0}, {}, {" <> StringRiffle[toks, ", "] <> "}, {" <>
            StringRiffle[netFacs, ", "] <> "}}", num, dr}
        ]] /@ opts;

(* Turn a component's factor list into the C++ Dirac-chain token(s) + its Lorentz "rest".
   Algorithm:
     1. pull out the Dirac heads (γ/γ5/σ/δ/dressed-numerator); if none, defer to compileTInv.
     2. order them into independent spinor loops (orderDiracLoops) and map each to a token string
        (tokenOf): a γ leg becomes dgamma/dslash, a σ becomes the matching dcomm* builder (legStr),
        a dressed-numerator becomes a SLOT (dtslot) whose DSlot is recorded.
     3. LOUD-GUARD that no Dirac head leaked into the non-Dirac "rest" (an un-collected dressed sum) —
        that would silently collapse/leak the trace; abort if so.
     4. compile the Lorentz rest separately (restCompiled = {restStr, scalar}) and return the Dirac
        core, scalar and rest SEPARATELY so the contraction stays deferred (rest emitted once, shared
        across a colour group). A dressed chain returns the ntDressedCore[chain, slots] marker instead. *)

compileDirac[factors_, ids_, env_, mask_, nc_] := Module[
    {diracFacs, loops, loopStrs, loopStrsBare, nEmptyLoops, slashVecs = {}, tokenOf, restFacs, restCompiled, legStr, sigStr, slots = {}, slotN = 0, dressed, vecOf},
(* μ -> the momentum q of an ntVec[q,μ] factor, built ONCE per call: tokenOf would otherwise rescan the
   whole factor list for every gamma token, which is quadratic in the factor count. Reverse before
   building the Association so a duplicate μ keeps the FIRST match, matching the old `First[...]`. *)
    vecOf = Association[Reverse[Cases[factors, ntVec[q_, m_] :> (m -> q)]]];
    diracFacs = Select[factors, MatchQ[#, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac | _ntDressedNum | _ntDiracSlot]&];
    If[diracFacs === {},
      Return[Append[compileTInv[Times @@ factors, ids, env, mask, nc], ""]]];
    loops = orderDiracLoops[diracFacs];(* one ordered token list per independent spinor loop *)
    dressed = !FreeQ[diracFacs, _ntDressedNum | _ntDiracSlot];
(* an ntSigma leg → its C++ arg string. A slashed leg is a list of {coeff, q} pairs (q a literal
   ntVec momentum, hence an env key) ⇒ a vlc `{{c1,b1},{c2,b2},…}` of (coeff, env Base) pairs (a
   single atomic momentum reduces to `{{1.,Base}}`); a free leg is its open Lorentz id. Guard loudly
   if a leg momentum is somehow not in env (it always should be — momentumOf collected it). *)
    legStr[{"slash", pairs_List}] := "{" <>
        StringRiffle[
          Function[cm,
              Module[{c = cm[[1]], q = cm[[2]]},
                If[!KeyExistsQ[env, q],
                  Print["[NumTracer] ERROR: σ slash leg momentum ", q, " absent from env ", Keys[env]];
                  Abort[]];
                "{" <> cppNum[c] <> "," <> ToString[env[q]["Base"]] <> "}"]
            ] /@ pairs,
          ", "
        ] <> "}";
    legStr[{"free", mu_}] := ToString[ids[mu]];
    (* the bare-commutator [A,B] builder, picked by each leg's free/slash kind (see et/inv/dirac.hpp). *)
    sigStr[legA_, legB_] := Module[{ta = legA[[1]], tb = legB[[1]], sa = legStr[legA], sb = legStr[legB]},
        Which[
          ta === "free" && tb === "free",
            "dcomm(" <> sa <> ", " <> sb <> ")",
          ta === "slash" && tb === "slash",
            "dcomm_ss(" <> sa <> ", " <> sb <> ")",
          ta === "free" && tb === "slash",
            "dcomm_fs(" <> sa <> ", " <> sb <> ")",
          True,
            "dcomm_sf(" <> sa <> ", " <> sb <> ")"]];
(* one chain token. A ntDressedNum becomes a SLOT (dtslot index, its DSlot recorded); the others are
   the usual gamma/slash/γ5/σ tokens. In the dressed case every fixed token is wrapped in dtfix(...)
   so the chain is a std::vector<DChainTok>; otherwise the bare token goes into a DiracNet{...}. *)
    tokenOf =
      Function[gf,
        Which[
          MatchQ[gf, _ntDressedNum],
            (
              AppendTo[slots, dressedSlotStr[gf, env]];
              With[{k = slotN},
                slotN++;
                "dtslot(" <> ToString[k] <> ")"]),
          MatchQ[gf, _ntDiracSlot],
            (
              AppendTo[slots, diracSlotStr[gf, ids, env, mask, nc]];
              With[{k = slotN},
                slotN++;
                "dtslot(" <> ToString[k] <> ")"]),
          MatchQ[gf, _ntGamma5],
            If[dressed,
              "dtfix(dg5())",
              "dg5()"],
          MatchQ[gf, _ntSigma],
            With[{s = sigStr[gf[[1]], gf[[2]]]},
              If[dressed,
                "dtfix(" <> s <> ")",
                s]],
          True,
            Module[{mu = First[gf], vq, t},
              vq = Lookup[vecOf, mu, Missing[]];
              t =
                If[!MissingQ[vq],
                  (
                    AppendTo[slashVecs, ntVec[vq, mu]];
                    "dslash({{1.0," <> ToString[env[vq]["Base"]] <> "}})"),
                  "dgamma(" <> ToString[ids[mu]] <> ")"];
              If[dressed,
                "dtfix(" <> t <> ")",
                t]]]];
    (* one token-string per spinor loop (mapped in order so the slot/slashVec side effects accumulate) *)
    loopStrs = (StringRiffle[tokenOf /@ #, ", "])& /@ loops;
(* COLLAPSED (token-free) spinor loops. orderDiracFacs drops δ "connector" factors, so a loop built
   ONLY from ntDeltaDirac (a bare closed spinor δ-loop, e.g. <P_2,T_2> of AqbqDirect8) yields an
   EMPTY token list. Downstream, ndetail::split_loops discards empty segments, so such a loop
   contributes NO factor at all and its tr(1) = 4 is silently lost — a 4x undercount per collapsed
   loop. (The DRESSED path already compensates for exactly this, restoring 4^(#collapsed loops); see
   the comment on nloops in numeric_value_dressed_netval. The bare path had no such compensation.)
   Fix it HERE rather than in C++: only the front end can tell a loop that collapsed to the identity
   (tr(1) = 4) from a component with no spinor loop at all (a pure-gauge diagram, where 4 must NOT be
   applied) — by the time the chain reaches the runtime both look like an empty DiracNet. compileDirac
   is only ever entered with diracFacs =!= {}, so every loop counted here is a genuine spinor loop.
   Flows whose every loop keeps a fixed γ never collapse => nEmpty == 0 => kernels byte-identical. *)
    nEmptyLoops = Count[loopStrs, ""];
(* REACHABILITY (checked 2026-07-18, and why the riffle below filters anyway): a token-free loop
   consists only of ntDeltaDirac, which carries NO Lorentz and NO colour label, so it can never
   share a label with another spinor loop — connectedComponents always isolates it into its own
   component. Hence loopStrs is exactly {""} whenever it holds an empty entry, and the empty entry
   is never adjacent to a non-empty one. Verified: a net with an all-δ loop, a γ loop and a colour
   loop splits into 3 components, the δ one reporting loops = {{}}.
   The filter is still applied because StringRiffle over a list MIXING "" with real tokens would
   emit stray commas — `DiracNet{, dloopsep(), dgamma(0)}` — which is a C++ syntax error in leading
   or middle position but a LEGAL trailing comma in last position. That is a correctness that
   depends on loop ORDER, which nothing enforces; it would become reachable the moment a Dirac head
   that carries a Lorentz index can also be token-free. Filtering makes the emission order-blind.
   Applied ONLY to the bare branch below: the DRESSED chain must keep every segment, because the
   runtime derives nCollapsed from the LoopSep marker count (ndetail::dress_enumerate) and
   dropping a segment would silently change it. A dressed loop is never token-free at codegen time
   anyway — an ntDressedNum is itself a chain token (a dtslot). *)
    loopStrsBare = DeleteCases[loopStrs, ""];
    restFacs = DeleteCases[factors, Alternatives @@ Join[diracFacs, slashVecs]];
(* LOUD GUARD against the silent dressed-numerator fall-through. `diracFacs` captures only BARE Dirac
   heads; a dressed propagator-numerator sum (Mq·δ − I·Z·γ·p̸) that was NEITHER distributed
   (expandBridges) NOR collected (rewriteDressedNums → ntDressedNum) survives as a `Plus`/`Power`
   factor with Dirac heads NESTED inside it. Such a factor lands in `restFacs` → `compileTInv`, which
   handles only Lorentz/colour/scalars: it then SILENTLY drops the nested γ/δ (the trace collapses,
   e.g. a pion loop → tr(γ5)=0) or leaks the raw head into the generated C++ (undeclared `ntGamma` /
   `ntDeltaDirac`). Either way the diagram is wrong with no warning (observed in the hSigL meson
   sector, ~1% of FORM). Detect it and abort with the offending factor rather than emit a wrong kernel.
   A correctly-handled diagram has all Dirac structure in `diracFacs` (bare heads or ntDressedNum), so
   `restFacs` is Dirac-free; this never trips on the validated flows. *)
    With[{leak = Select[restFacs, !FreeQ[#, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac | _ntDressedNum | _ntDiracSlot]&]},
      If[leak =!= {},
        Print["[NumTracer] ERROR: un-handled Dirac structure in a non-Dirac factor — a dressed ", "propagator-numerator sum was NEITHER distributed NOR collected into ntDressedNum, so the ", "numeric backend would silently drop/leak its gamma structure (collapsed trace or ", "untranslated C++). This is a front-end collection gap (collectibleDiracSumQ rejected a sum ", "that distributeQ also skipped). Offending factor(s):\n  ", leak];
        Abort[]]];
    restCompiled = compileTInv[Times @@ restFacs, ids, env, mask, nc];(* {restStr, scal}; restStr "" if trivial *)
(* return the Dirac core and the Lorentz "rest" (projector) SEPARATELY: {coreStr, scal, restStr}.
   The contract is deferred — splitColourGroupsInv factors a rest shared across a colour group out of
   the Dirac-trace sum so it is emitted ONCE (a shared sub-net) and contracted lazily in the
   generator, never materialising the |trace|×|projector| product. restStr=="" ⇒
   no projector. The numeric backend emits the bare DiracNet{...}; the runtime contracts it against
   the Lorentz rest via numeric_value_netval (the matrix-product trace). A DRESSED chain (one or more
   ntDressedNum slots) is returned as the marker ntDressedCore[chainStr, slotsStr]; the generator
   contracts it via numeric_value_dressed_netval_mp — a plain STRUCTURAL MPoly, since lever (b) strips
   each option's dressing into a per-sub-term scalar + monomial before contraction, and the DPoly is
   assembled later in phase B (fold_groups_streaming_dressed). *)
(* join the per-loop token strings with a loop separator so the runtime traces each spinor loop
   independently and contracts their shared gluon legs; a single-loop component has no separator and
   emits the identical net as before. *)
    If[dressed,
(* Carry the slots STRUCTURED (a list over chain slots of that slot's option-string list) rather than
   pre-joined into one multi-option "std::vector<DSlot>{…}" string. emitNumericGenerator expands the
   Cartesian product of the options into one single-option dressed sub-term per combination, so each
   combination becomes an ordinary trace that dedups across nets and contracts over phase A's flat
   parallel work list — instead of C++ dress_collect enumerating the 3ⁿ combinations serially per net. *)
      Module[{chain = "std::vector<DChainTok>{" <> StringRiffle[loopStrs, ", dtfix(dloopsep()), "] <> "}"},
        {ntDressedCore[chain, slots], restCompiled[[2]], restCompiled[[1]]}],
      Module[{core = "DiracNet{" <> StringRiffle[loopStrsBare, ", dloopsep(), "] <> "}"},
        {core, restCompiled[[2]] * 4^nEmptyLoops, restCompiled[[1]]}]]];

(* ---- numeric (matrix-product) backend: component table + user symbols (task #22) -------------
   The numeric backend has NO sp-invariant basis: it evaluates scalar products numerically from each
   momentum's 4 COMPONENTS. So instead of an sp-invariant decomposition it needs
   only: the polynomial variables (the free user symbols), each fundamental momentum's 4 components as
   MPoly-builder C++, and the C++ fill formula for each symbol (a kernel arg, or a derived symbol like
   sin1 = sqrt(1-cos1^2)). Composite momenta resolve by component arithmetic via resolveComponents,
   exactly like the frame path — the user "Components" association is just a frame whose entries are
   polynomial. *)
(* Polynomialise a frame whose loop direction carries Sqrt[1-cos^2]: introduce a sin symbol per angle
   and record its definition. A user-supplied (already polynomial) "Components" assoc is a no-op. *)
(* Unit-loop spec (the COMPACT numeric parametrisation that matches the sp invariant count). The
   naive frame fallback bloats because the LOOP momentum is decomposed into angle PRODUCTS
   (l1·cos1, l1·sin1·cos2, … — degree-3), so every l·p_i expands into many angle monomials and the
   degrees compound through a projector chain. Instead write the loop as MAGNITUDE × UNIT-DIRECTION:
   comp_μ = l1 · Uμ where Uμ are degree-1 symbols (fill computes Uμ = dirμ from the kernel angles).
   Then l·p_i = l1·p·Σ Uμ u_i^μ is degree-1 in Uμ (as compact as sp(l,p_i)), and the bare-loop
   denominator Σ(l1·Uμ)² = l1²·ΣUμ² collapses to the MONOMIAL l1² under the unit constraint ΣUμ²=1
   (returned as a unit `group`), so its 1/l² atom cancels — exactly like inv's `rel`. Externals depend
   only on `p` (kept as the numeric p-vector). The loop is the momentum whose components carry `magSym`. *)

unitLoopFrameSpec[frame_, pSym_, magSym_] := Module[{defs = <||>, groups = {}, n = 0, nf},
    nf =
      Association @
        KeyValueMap[
          Function[{q, comps},
            q ->
              If[SubsetQ[{pSym}, Variables[comps]],
                comps,
                (* external: numeric p-vector *)
                Module[
                  {grp = {}, nc},
                  (* loop: comp_μ = magSym · Uμ *)
                  nc =
                    Table[
                      Module[{dir = Coefficient[comps[[mu]], magSym], s},
                        If[dir === 0,
                          0,
                          s = Symbol["ntU$" <> ToString[n++]];
                          defs[s] = dir;
                          AppendTo[grp, s];
                          magSym s]],
                      {mu, 1, 4}];
                  AppendTo[groups, grp];
                  nc]]],
          frame];
    {nf, defs, groups}];

(* General frame -> {polyFrame, defs}: replace the trig sub-expressions in the frame components with
   fresh symbols so the components become polynomials in (magnitudes, cos's, sin's). Introduces a
   symbol for each loop polar factor Sqrt[1-cos^2] (a sin) and each bare-angle Cos/Sin (the 4-point
   frame's φ); `defs` maps each new symbol back to its closed form (so the kernel can fill it). The
   fallback path when the unit-loop parametrisation (unitLoopFrameSpec) does not apply. *)

polyFrameSpec[frame_] := Module[{defs = <||>, vals, rules = {}, polyFrame},
    vals = PowerExpand[Values[frame]];
    (* loop polar factor Sqrt[1-cos^2] -> a sin symbol *)
    Do[
      With[{s = Symbol["ntSin$" <> SymbolName[c]]},
        defs[s] = Sqrt[1 - c^2];
        AppendTo[rules, Sqrt[1 - c^2] -> s]],
      {c, DeleteDuplicates @ Cases[vals, Power[1 - c_^2, 1/2] :> c, Infinity]}];
    (* trig of a bare angle symbol (e.g. Cos[phi], Sin[phi] in the 4-point frame) -> cos/sin symbols *)
    Do[
      With[{s = Symbol["ntCos$" <> SymbolName[a]]},
        defs[s] = Cos[a];
        AppendTo[rules, Cos[a] -> s]],
      {a, DeleteDuplicates @ Cases[vals, Cos[a_Symbol] :> a, Infinity]}];
    Do[
      With[{s = Symbol["ntSinA$" <> SymbolName[a]]},
        defs[s] = Sin[a];
        AppendTo[rules, Sin[a] -> s]],
      {a, DeleteDuplicates @ Cases[vals, Sin[a_Symbol] :> a, Infinity]}];
    polyFrame = Association @ Thread[Keys[frame] -> (vals /. rules)];
    {polyFrame, defs}];

(* Whether `frame` matches the unit-loop spec's assumption: every momentum is either an external
   depending only on `pSym`, or a loop whose every component is `dir·magSym` (proportional to the
   single magnitude). A finite-T frame breaks this — the external carries an independent temporal
   component p0 and the loop an independent l0 (neither ∝ magSym) — so we must fall back to the
   general polyFrameSpec. *)

unitLoopOkQ[frame_, pSym_, magSym_] := AllTrue[
    Values[frame],
    Function[comps,
      Module[{cc = PowerExpand[comps]},
        SubsetQ[{pSym}, Variables[cc]] || AllTrue[Range[4], (Simplify[cc[[#]] - Coefficient[cc[[#]], magSym] magSym] === 0)&]
      ]]];

numericComponents[env_, frame_, symDefs_, unitGroups_ : {}] := Module[
    {compExpr, usyms, nsym, mpcpp, compCpp, varFill, vfill, idx, units},
    (* 4 components per momentum Base (polynomial in the user symbols). *)
    compExpr = Association @ KeyValueMap[#2["Base"] -> resolveComponents[#1, frame]&, env];
(* reject only fractional powers of NON-numeric bases (e.g. Sqrt[1-cos1^2] not rewritten to a sin
   symbol); a numeric irrational coefficient like Sqrt[3]/2 from the 120-degree external frame is a
   perfectly good polynomial coefficient. *)
    With[{bad = DeleteDuplicates @ Cases[Values[compExpr], Power[b_, _Rational] /; !NumericQ[b], Infinity]},
      If[bad =!= {},
        Print["numericComponents: non-polynomial components (fractional power of a symbol remains): ", bad];
        Abort[]]];
    usyms = Sort @ DeleteDuplicates @ Flatten[Variables /@ Values[compExpr]];
    nsym = Length[usyms];
(* Emit through the generator's `env` (a LorentzEnv bound to nsym) — the sole construction path for
   MPoly now that the bare-nsym factories are private. *)
    mpcpp[e_] := Module[{rules = CoefficientRules[e, usyms]},
        If[rules === {},
          "env.zero()",
          "(" <> StringRiffle[("env.mono({" <> StringRiffle[ToString /@ #[[1]], ","] <> "},Cx{" <> cppNum[#[[2]]] <> ",0})")& /@ rules, " + "] <> ")"
        ]];
    compCpp = Association @ KeyValueMap[#1 -> (mpcpp /@ #2)&, compExpr];
    vfill[s_] := If[KeyExistsQ[symDefs, s],
        cppFlat[symDefs[s]],
        SymbolName[s]];
    varFill = vfill /@ usyms;(* indexed by MPoly var id (0-based) *)
(* Unit-constraint groups (ΣUμ²=1) as MPoly var-index lists, so the C++ reduce_units collapses the
   bare-loop denominator to the monomial l1² and the U·U projector factors. The caller supplies the
   groups as symbol lists (the loop's unit-direction components); drop any symbol not in usyms (a
   component that vanished) and any group with <2 surviving symbols. *)
    idx = AssociationThread[usyms -> Range[Length[usyms]] - 1];
    units =
      DeleteCases[
        Function[g,
            Lookup[idx, Select[g, KeyExistsQ[idx, #]&], Nothing]
          ] /@ unitGroups,
        _ ? (Length[#] < 2&)];
    <|
      "nsym" -> nsym,
      "usyms" -> usyms,
      "compCpp" -> compCpp,
      "varFill" -> varFill,
(* default -1 (-> comp size 0) for a purely scalar integrand whose 4-vector component env is
   empty: Max[{}] is -Infinity, which would leak into the C++ as the comp() vector size. The loop
   momentum of such a flow (e.g. a bosonic meson-potential tadpole) is carried by usyms/units. *)
      "symNamesCpp" -> (("(" <> # <> ")")& /@ varFill),
      "maxBase" -> Max[Append[#["Base"]& /@ Values[env], -1]],
      "units" -> units
    |>];

(* ---- emit the NUMERIC generator program. Builds each net's DiracNet chain + Lorentz/projector
        NetVal in C++, contracts them numerically (4×4 matrix products), folds colour per group, and
        PRINTS the committed straight-line kernel header. No reduce/rebase/ibp/sp-kinematics. *)

emitNumericGenerator[invNets_, invRest_, colourNets_, groups_, ncomp_, nsInner_, fillArgSig_, kns_:"numtracer_kernels", complexQ_:False, realOnlyG_ : {}, crossCSE_:False] :=
  Module[{nNet = Length[invNets], nGrp = Length[groups], nsym = ncomp["nsym"], maxBase = ncomp["maxBase"], varFill = ncomp["varFill"], symNames = ncomp["symNamesCpp"], compCpp = ncomp["compCpp"], unitG = ncomp["units"], str, tmpl, pre, unitPre, nUnits, units, decl, diracNetStrs, lorentzNetStrs, subScalars, dressChains, dressSlotOpts, hasDressed, allDefs, main, compInit, cseDefs, cseDecls, chunkDecls, ntNoDedup, subKeys, netTerms, refCount, distinctSubs, subIdxOf, nSub, nReused, sdnDefs, slnDefs, sdchDefs, sdslDefs, sdnCDecl, slnCDecl, sdchCDecl, sdslCDecl,
    chpDefs, chpCDecl, optpDefs, optpCDecl, sdchrDefs, sdchrCDecl, sdslrDefs, sdslrCDecl, dressAtomIds},
    str[x_] := ToString[x];
(* DRESSED nets (symbolic dressing collection): a core may be ntDressedCore[chainStr, slotsStr]
   (a numerator structure-sum kept eager). LEVER (b): the trace table is PLAIN MPoly for both paths — a
   dressed sub-term's structural trace (dressing stripped) contracts via numeric_value_dressed_netval_mp,
   and its dressing rides the per-sub-term scalar (dsc, numeric) + monomial (sdr, atom ids). The phase-B
   fold (fold_groups_streaming_dressed) assembles those into the per-net DPoly, and the dressings ride the
   env as kind-2 `dress` leaves filled by fm.dress. So combinations that share a concrete structure but
   differ only in dressing collapse to ONE trace; the non-dressed path is byte-identical. *)
    hasDressed = !FreeQ[invNets, _ntDressedCore];
    (* shared wrapper templates (same as the inv generator) so the net-builder strings compile. *)
    tmpl = "template<int Mu,int Nu,int Lb,int Mask,int Inv> NetVal tproj(){ return projT(Mu,Nu,Lb,Inv); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int Inv> NetVal lproj(){ return projL(Mu,Nu,Lb,Inv); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int InvS> NetVal mproj(){ return projM(Mu,Nu,Lb,InvS); }\n" <> "template<int Mu,int Nu,int Lb,int Mask,int Inv,int InvS> NetVal eproj(){ return projE(Mu,Nu,Lb,Inv,InvS); }\n" <> "template<int Mu,int Nu> NetVal lmetric(){ return met(Mu,Nu); }\n" <> "template<int Lbl,int Base,int Mask> NetVal lvec(){ return vec(Lbl,Base); }\n" <> "template<int A,int B,int C,int D> NetVal leps(){ return epsilon(A,B,C,D); }\n" <> "inline NetVal konst(double c){ return NetVal{PTerm{Cx{c,0}, {}}}; }\n" <> "template<class L> struct litco;\n" <> "template<numtracer::Cx C> struct litco<numtracer::Lit<C>>{ static constexpr numtracer::Cx v=C; };\n" <> "template<class L> NetVal sc(NetVal x){ return scale(litco<L>::v, std::move(x)); }\n";
(* per net: a colour group is a SUM of sub-terms. invNets[i] = {core_b…} (each a DiracNet literal
   for a gamma branch, or a Lorentz NetVal for a gamma-free branch); invRest[i] = {{rest_b,scal_b}…}
   parallel. Build per net the parallel lists of {DiracNet builder, NetVal builder} (a gamma-free
   branch → empty DiracNet + the whole net as the rest) plus the sub-term scalars; the generator
   sums mp[i] = Σ_b scal_b · numeric_value_netval(dn[i][b], ln[i][b]). *)
(* Each branch yields a LIST of {ds, ls, scal, dc, dl, dr} sub-terms (usually length 1). A DRESSED branch
   expands the Cartesian product of its chain's slot options (Tuples) into ONE structural sub-term per
   combination. LEVER (b): each slot option is now a TRIPLE {structStr, num, dr} (dressing-free DSlotOpt +
   numeric Cx + dress-atom ids). For each combination we keep the STRUCTURE (dl = the list of structStr,
   one per slot) for the trace key, fold the numeric part (∏ num) INTO the sub-term scalar, and carry the
   dressing-atom multiset (dr = ⋃ dress) as the DPoly key. So combinations that share a concrete structure
   but differ only in dressing collapse to ONE plain-MPoly trace (the 6.2× dedup), and the trace table
   loses its dressing dimension entirely. A slot's option list is nv[[2]][[k]]; Tuples over them gives every
   combination. The non-dressed branches carry an empty dress key, so non-dressed flows stay byte-identical. *)
    {diracNetStrs, lorentzNetStrs, subScalars, dressChains, dressSlotOpts, dressAtomIds} =
      Transpose @
        MapThread[
          Function[{cores, rss},
            If[cores === {},
              {{}, {}, {}, {}, {}, {}},
              Transpose[
                Join @@
                  MapThread[
                    Function[{nv, rv, scal},
                      Module[{lsStr = If[rv === "", "NetVal{}", rv]},
                        Which[
                          MatchQ[nv, _ntDressedCore],(* dressed numerator: expand slot options → structural sub-terms *)
                            With[{chain = nv[[1]], slotOpts = nv[[2]]},
                              (* dl (5th) = this combination's STRUCTURAL option-string LIST (one dressing-free
                                                structStr per chain slot), kept STRUCTURED so the table emitter pools the distinct
                                                structures; the numeric Cx folds into scal and the dress ids become dr (6th). *)
                              Function[combo,
                                {"DiracNet{}", lsStr,
                                scal * (Times @@ (#[[2]]& /@ combo)),   (* fold ∏ numeric Cx into the scalar *)
                                chain,
                                #[[1]]& /@ combo,                        (* structural options (dressing-free) *)
                                Sort[Catenate[#[[3]]& /@ combo]]}        (* ⋃ dress atom ids, sorted (Catenate[{}]={}) *)
                              ] /@ If[slotOpts === {}, {{}}, Tuples[slotOpts]]],
                          StringMatchQ[nv, "DiracNet" ~~ ___],(* gamma branch: DiracNet + projector rest *)
                            {{nv, lsStr, scal, "std::vector<DChainTok>{}", {}, {}}},
                          True,(* gamma-free branch: whole net is the rest *)
                            {{"DiracNet{}", nv, scal, "std::vector<DChainTok>{}", {}, {}}}]]],
                    {cores, rss[[All, 1]], rss[[All, 2]]}]]]],
          {invNets, invRest}];
    pre =
      StringJoin[
        "// GENERATED by MakeNTKernel — do not edit. Numeric matrix-product tensor traces.\n",
(* one umbrella header pulls the whole engine API (network/dirac/gen/sun_net/numeric_*/core) for
   this single generator TU. The parallel -O0 net-builder units below deliberately keep their
   minimal includes — pulling the umbrella into each would re-parse the whole engine per unit and
   regress generation time. The emitted kernel stays minimal too (runtime.hpp + sun_data.hpp). *)
        "#include \"numtracer/numtracer.hpp\"\n",
(* the two contraction phases (phase A: contract each distinct trace once over a flat work list;
   phase B: fold each net's traces as a balanced tree). Host-only (it spawns threads) and MAIN-TU
   ONLY — the net-builder units are compiled -fno-exceptions and must not see it. *)
        "#include \"numtracer/numeric/trace_fold.hpp\"\n",
        "#include <iostream>\n#include <string>\n#include <utility>\n#include <vector>\n#include <array>\n#include <cstdlib>\n",
        "#include <thread>\n#include <atomic>\n#include <mutex>\n#include <chrono>\n#include <cstdio>\n#include <algorithm>\n#include <system_error>\n",
        "#include <unistd.h>\n#if defined(__GLIBC__)\n#include <malloc.h>\n#endif\n",
(* Resident set in MB, straight from /proc — so the generator can report its OWN peak rather than
   the caller having to wrap it in `/usr/bin/time`. Returns 0 where /proc is absent. *)
        "static double ntRssMB(){ long pages=0; if(FILE* f=std::fopen(\"/proc/self/statm\",\"r\")){ long tot=0; if(std::fscanf(f,\"%ld %ld\",&tot,&pages)!=2) pages=0; std::fclose(f); } return pages*(double)sysconf(_SC_PAGESIZE)/1048576.0; }\n",
        "#include <sstream>\n#include <unordered_map>\n",
        "using numtracer::Cx;\n",
        "namespace numtracer::network {\n",
        tmpl,
        "}\n",
        "using namespace numtracer::network;\nusing namespace numtracer::numeric;\n"];
    unitPre =
      "// GENERATED by MakeNTKernel — do not edit. Numeric net-builder unit (compiled -O0).\n" <> "#include \"numtracer/network/network.hpp\"\n#include \"numtracer/network/dirac.hpp\"\n#include \"numtracer/core/lit.hpp\"\n#include <utility>\n" <>
(* dressed nets emit dch<i>()/dsl<i>() builders here (the big DChainTok/DSlot literals — moved OFF
   the single -O1 main TU onto these parallel -O0 units, since a 100k+-char braced-init is ~quadratic
   even at -O0); they need the dressed-token types from numeric_contract.hpp. Non-dressed units don't
   include it (stay byte-identical + fast). *)
        If[hasDressed,
          "#include \"numtracer/numeric/numeric_contract.hpp\"\n",
          ""
        ] <> "using numtracer::Cx;\nnamespace numtracer::network {\n" <> tmpl <> "}\nusing namespace numtracer::network;\n" <>
        If[hasDressed,
          "using namespace numtracer::numeric;\n",
          ""];
(* NET-LEVEL CSE: dense projections (e.g. the σ^μν struct-7 quark-gluon vertex) emit the SAME net
   sub-term thousands of times — measured 48624 lnet terms but only 66 distinct (534x), so the raw
   generator C++ balloons to ~25 MB and the -O0 compile dominates generation. Hash-cons each DISTINCT
   net term into a shared accessor `lc<k>()` / `dc<k>()` (a function-local `static const` so the net
   is also BUILT once at run time, not once per use), and reference it. Trivial/empty literals and
   unique terms stay inline. Correctness-preserving: each use copies the shared NetVal/DiracNet,
   exactly as the inlined expression did. *)
    cseDefs = {};
    cseDecls = "";
(* The lc<k>()/dc<k>() net-level CSE runs for BOTH paths: the dressed (collected) lnet builders are
   just as repetitive as the dense σ case (the same Lorentz/projector structures), and the accessors
   only change the GENERATOR's internal sharing — each use copies the shared NetVal/DiracNet, so the
   contracted result (hence the committed kernel) is byte-identical. Previously skipped when dressed,
   which left the dressed generator C++ bloated and its -O0 compile ~2x slower. *)
    Module[{lTerms = Flatten[lorentzNetStrs], dTerms = Flatten[diracNetStrs], lCount, dCount, lMap = <||>, dMap = <||>, li = 0, di = 0, lRef, dRef},
      With[{
        ntT =
          First @
            AbsoluteTiming[
              lCount = Counts[lTerms];
              dCount = Counts[dTerms];]},
        ntLog["[prof] CSE Counts (", Length[lTerms], "+", Length[dTerms], " terms): ", ntT, " s"]];
      KeyValueMap[
        Function[{t, c},
          If[c >= 2 && t =!= "NetVal{}" && t =!= "",
            lMap[t] = "lc" <> str[li];
            li++]],
        lCount];
      KeyValueMap[
        Function[{t, c},
          If[c >= 2 && t =!= "DiracNet{}" && t =!= "",
            dMap[t] = "dc" <> str[di];
            di++]],
        dCount];
      lRef[t_] := If[KeyExistsQ[lMap, t],
          lMap[t] <> "()",
          t];
      dRef[t_] := If[KeyExistsQ[dMap, t],
          dMap[t] <> "()",
          t];
      cseDefs =
        Join[
          KeyValueMap[
            Function[{t, nm},
              "const DiracNet& " <> nm <> "(){ static const DiracNet v = " <> t <> "; return v; }"],
            dMap],
          KeyValueMap[
            Function[{t, nm},
              "const NetVal& " <> nm <> "(){ static const NetVal v = " <> t <> "; return v; }"],
            lMap]];
      cseDecls =
        StringJoin[
          Riffle[
            Join[
              KeyValueMap[
                Function[{t, nm},
                  "const DiracNet& " <> nm <> "();"],
                dMap],
              KeyValueMap[
                Function[{t, nm},
                  "const NetVal& " <> nm <> "();"],
                lMap]],
            "\n"]];
      With[{
        ntT =
          First @
            AbsoluteTiming[
              diracNetStrs = Map[dRef, diracNetStrs, {2}];
              lorentzNetStrs = Map[lRef, lorentzNetStrs, {2}];]},
        ntLog["[prof] CSE ref-rewrite: ", ntT, " s"]];
      ntLog["[cse] net terms: lnet ", Length[lTerms], "->", Length[lMap], " distinct, dnet ", Length[dTerms], "->", Length[dMap], " distinct shared builders"]
    ];
(* ---- GLOBAL SUB-TERM DEDUP ------------------------------------------------------------------
   A net is Σ_b scal_b · contract(dn_b, ln_b, dch_b, dsl_b). The generator's cost is one trace
   contraction per (net, sub-term) — but the SAME (dn,ln,dch,dsl) tuple recurs across nets and
   colour branches, so most of those contractions recompute a trace already computed. Measured on
   the dense flows: 30,807 contractions for 6,041 distinct traces (5.1x), and 246,456 for 32,784
   (7.5x). So contract each DISTINCT trace ONCE into a shared table and let every net fold the
   table with its own scalars. Two independent wins:
     - the contraction phase (the bottleneck) shrinks by the redundancy factor;
     - the parallel phase becomes a FLAT list of uniform work items. Scheduling per NET could not
       use the machine: sub-terms per net are wildly skewed (max 2880 vs a median of 27), so the
       single biggest net alone exceeded the ideal per-thread load and pinned utilisation at ~33%
       however many cores were available.
   Sub-terms sharing a trace are merged and their scalars SUMMED (which also shortens each net's
   fold), and a merged term whose scalars sum to 0 is dropped (σ-commutator cancellations).

   This SUPERSEDES the old per-net (dn,ln) merge, which ran only when !hasDressed on the theory that
   "the collected path emits ~1 sub-term per net, so there is nothing to merge". That was false for
   dense flows (187 sub-terms/net). The global key reduces to (dn,ln) on the non-dressed path (where
   dch/dsl are constant ""), so it does everything the old merge did, and across nets as well.

   NT_GEN_NO_DEDUP=1 turns the dedup off: every occurrence becomes its own trace, nothing is merged
   or dropped, and nothing is cached (nReused=0), so each net contracts its own sub-terms on demand
   — the pre-dedup behaviour. It is the escape hatch, and the control for the equivalence test
   (generate twice, compare the kernels' VALUES; a byte-diff is meaningless because dedup changes
   GlobalEnv interning order and so renumbers every sN). *)
    ntNoDedup = (Environment["NT_GEN_NO_DEDUP"] =!= $Failed);
(* LEVER (b): each sub-term carries a {traceKey, dressKey} PAIR. The traceKey is dressing-free
   {ds, ls, chain, structural-options} — so combinations that share a concrete structure but differ only
   in dressing collapse to ONE trace. The dressKey is the sorted dress-atom multiset; it becomes the
   DPoly channel a net's fold routes this sub-term's (scaled) trace into. Sub-terms MERGE only when they
   share BOTH the trace AND the dress channel (a merge across channels would corrupt the DPoly). For a
   non-dressed flow every dressKey is {}, so the grouping reduces to the old traceKey grouping and the
   emitted kernel is byte-identical. *)
    subKeys =
      MapThread[
        Function[{ds, ls, dc, dl, dr, ni},
          Table[
            {
              Join[
                {
                  ds[[j]],
                  ls[[j]],
                  If[hasDressed,
                    dc[[j]],
                    ""],
                  If[hasDressed,
                    dl[[j]],
                    ""]},
                If[ntNoDedup,
                  {ni, j},
                  {}]],
              (* dressKey: the dress-atom multiset ({} for a non-dressed sub-term) *)
              dr[[j]]},
            (* unique per occurrence => nothing ever merges *)
            {j, Length[ds]}]],
        {
          diracNetStrs,
          lorentzNetStrs,
          If[hasDressed,
            dressChains,
            diracNetStrs],
          If[hasDressed,
            dressSlotOpts,
            diracNetStrs],
          dressAtomIds,
          Range[Length[diracNetStrs]]}];
(* Per net: merge the sub-terms that share a trace AND a dress channel, summing scalars; drop the zero
   sums. Do this BEFORE counting references — a (trace, channel) occurring twice inside ONE net collapses
   to a single reference here, so its raw occurrence count would overstate its reuse, and a channel whose
   scalars cancel is not referenced at all and must not be contracted. *)
    netTerms =
      MapThread[
        Function[{keys, scs},
          (* GatherBy element = {{traceKey, dressKey}, sc}; g[[1,1]] is the {tk,dk} PAIR, so tk = g[[1,1,1]],
                  dk = g[[1,1,2]]. The trace table keys on tk ALONE (refCount below), so dressing variants of one
                  structure share the trace slot; only the fold entry carries dk. *)
          Select[(Function[g, {g[[1, 1, 1]], g[[1, 1, 2]], Total[g[[All, 2]]]}]) /@ GatherBy[Transpose[{keys, scs}], First],
            #[[3]] =!= 0&]],
        {subKeys, subScalars}];
(* Order the distinct traces by DESCENDING reference count, so a memory-capped run still caches the
   traces that repay caching most, and the singletons (refCount 1 — computing one costs the same
   whether or not it is cached, so caching it is pure RAM for no saving) land at the end where the
   default cap excludes them. *)
    (* First /@ #, not #[[All,1]] — the latter errors on a net whose terms all cancelled to {} *)
    refCount = Counts[Flatten[Map[First /@ #&, netTerms], 1]];
    distinctSubs = Keys[ReverseSort[refCount]];
    subIdxOf = AssociationThread[distinctSubs -> Range[0, Length[distinctSubs] - 1]];
    nSub = Length[distinctSubs];
(* Published for the COMPILE step (mainOpt), which runs later, in another function, and needs the
   flow's distinct-trace count to pick the main-TU optimisation level — the generator RUN scales
   with it. A global rather than a threaded argument because the two are already strictly
   sequential within one generation; reset per generation right here, where it is computed, so a
   stale value from a previous flow can never leak into the next one's decision. *)
    $ntGenNSub = nSub;
    nReused = Count[Values[refCount], c_ /; c >= 2];(* == the length of the sorted prefix worth caching *)
    netTerms =
      Map[
        Function[nt,
            {subIdxOf[nt[[1]]], nt[[2]], nt[[3]]}
          ] /@ #&,
        netTerms];
    ntLog[
      "[cse] sub-terms: ",
      Total[Length /@ subKeys],
      " contractions -> ",
      nSub,
      " distinct traces (",
      ToString @ NumberForm[N[Total[Length /@ subKeys] / Max[1, nSub]], {5, 2}],
      "x), ",
      nReused,
      " reused (cached), ",
      nSub - nReused,
      " singletons; per-net folds total ",
      Total[Length /@ netTerms],
      " (longest ",
      Max[Append[Length /@ netTerms, 0]],
      ")",
      If[ntNoDedup,
        " [NT_GEN_NO_DEDUP]",
        ""]];
(* The distinct traces are emitted as ONE flat table each (still chunked by ntChunkDef, whose helpers
   the bin-packer scatters across the -O0 units); the nets reference them by index. Guard nSub == 0:
   {}[[All,1]] is an error, not an empty list. *)
    {sdnDefs, sdnCDecl} =
      ntChunkDefs[
        "sdn",
        "std::vector<DiracNet>",
        If[nSub === 0,
          {{}},
          {distinctSubs[[All, 1]]}]];
    {slnDefs, slnCDecl} =
      ntChunkDefs[
        "sln",
        "std::vector<NetVal>",
        If[nSub === 0,
          {{}},
          {distinctSubs[[All, 2]]}]];
(* DRESSED slot tables — INTERNED (chain pool + option pool + per-sub-term index arrays).
   Expanding each structure×dressing COMBINATION into its own single-option sub-term (so phase A
   contracts them in parallel, killing the serial dress_collect) makes the chain and slot columns
   MASSIVELY redundant: a net's thousands of combinations share ONE chain (za3_147: 5 distinct over
   12101 sub-terms) and draw their DSlotOpts from a tiny per-slot option pool (35 distinct over 84672
   emissions — 99.9% redundant). Emitting the full literals per sub-term blew the generator SOURCE to
   ~9.5 MB (compile 2.7 s -> 13.5 s). Instead emit the distinct chains (`chp`) and options (`optp`)
   ONCE, and each sub-term as compact INDICES (`sdchR`: its chain index; `sdslR`: one option index per
   slot); main rebuilds sdch[k]/sdsl[k] from them in an O(nSub) loop — the exact hash-consing the
   sidx/dsc tables already use. Source ~9.5 MB -> ~0.4 MB, so the phase-A run-time win no longer costs
   a compile-time regression. Non-dressed sub-terms carry an empty chain / empty option list (the
   sdch[k].empty() numeric_value_netval fast path is preserved). *)
    {chpDefs, chpCDecl, sdchrDefs, sdchrCDecl, optpDefs, optpCDecl, sdslrDefs, sdslrCDecl} =
      If[hasDressed,
        Module[{chainStrs, combos, uChains, chainPos, uOpts, optPos, chainDefs, chainDecl, chainRefDefs, chainRefDecl, optDefs, optDecl, slotRefDefs, slotRefDecl},
          chainStrs = If[nSub === 0, {}, distinctSubs[[All, 3]]];
          combos    = If[nSub === 0, {}, distinctSubs[[All, 4]]];(* per sub-term: its option-string LIST ({} for a non-slot sub-term) *)
          uChains = DeleteDuplicates[chainStrs];
          chainPos = AssociationThread[uChains -> Range[Length[uChains]] - 1];
          uOpts = DeleteDuplicates[Flatten[combos]];
          optPos = AssociationThread[uOpts -> Range[Length[uOpts]] - 1];
          {chainDefs, chainDecl} = ntChunkDefs["chp", "std::vector<std::vector<DChainTok>>", {uChains}];
          {chainRefDefs, chainRefDecl} = ntChunkDefs["sdchR", "std::vector<int>", {ToString /@ (chainPos /@ chainStrs)}];
          {optDefs, optDecl} = ntChunkDefs["optp", "std::vector<DSlotOpt>", {uOpts}];
          {slotRefDefs, slotRefDecl} = ntChunkDefs["sdslR", "std::vector<std::vector<int>>",
            {("{" <> StringRiffle[ToString /@ (optPos /@ #), ","] <> "}")& /@ combos}];
          {chainDefs, chainDecl, chainRefDefs, chainRefDecl, optDefs, optDecl, slotRefDefs, slotRefDecl}],
        {{}, "", {}, "", {}, "", {}, ""}];
    chunkDecls = sdnCDecl <> slnCDecl <> chpCDecl <> sdchrCDecl <> optpCDecl <> sdslrCDecl;
    allDefs = Join[cseDefs, sdnDefs, slnDefs, chpDefs, sdchrDefs, optpDefs, sdslrDefs];
(* Size-aware unit count (~$ntUnitChars per unit, 8..$ntUnitCap): the CSE accessors are many but
   small, so the old 12-defs/unit rule would emit hundreds of tiny TUs each re-parsing the shared
   decl header. Defs are packed into the units by GREEDY BIN-PACKING (largest def first, into the
   currently-smallest unit — LPT) rather than round-robin by index: def sizes are skewed, and
   round-robin left a large max/median spread that made the biggest units the makespan stragglers.
   The cap must stay above what the size target asks for, or units quietly grow back past it. Now that
   ntChunkDef bounds any SINGLE def to ~$ntDefChunk, bin-packing can actually hit the target: a single
   oversized def is no longer an irreducible floor. *)
(* Total[StringLength/@...], NOT StringLength[StringJoin[...]] — the latter materialises every def as
   one giant string just to measure it (much slower, and it allocates the lot). *)
    nUnits = Min[Min[$ntUnitCap, Max[8, Ceiling[Total[StringLength /@ allDefs] / $ntUnitChars]]], Max[1, Length[allDefs]]];
    units =
      If[allDefs === {},
        {},
        (unitPre <> StringRiffle[#, "\n"] <> "\n")& /@
          Module[{bins = ConstantArray[{}, nUnits], loads = ConstantArray[0, nUnits], b},
            Do[
              b = First @ Ordering[loads, 1];
              bins[[b]] = Append[bins[[b]], d];
              loads[[b]] += StringLength[d],
              {d, allDefs[[Reverse @ Ordering[StringLength /@ allDefs]]]}];
            bins]];
(* the shared decl header: net-builder + CSE-accessor forward declarations (the units that call the
   lc<k>()/dc<k>() accessors #include this — emitted ONCE here, not duplicated per unit). *)
    decl =
      "// GENERATED by MakeNTKernel — do not edit. Numeric net-builder declarations.\n#pragma once\n" <> "#include \"numtracer/network/network.hpp\"\n#include \"numtracer/network/dirac.hpp\"\n#include <vector>\n" <>
        If[hasDressed,
          "#include \"numtracer/numeric/numeric_contract.hpp\"\n",
          ""
        ] <> "using namespace numtracer::network;\n" <>
        If[hasDressed,
          "using namespace numtracer::numeric;\n",
          ""
        ] <> cseDecls <> "\n" <>
(* the DISTINCT-trace tables (see the global sub-term dedup above): one flat builder each, which
   the nets index into — not one builder per net, as before the dedup. *)
        "std::vector<DiracNet> sdn0();\n" <> "std::vector<NetVal> sln0();\n" <>
        If[hasDressed,
          "std::vector<std::vector<DChainTok>> chp0();\n" <> "std::vector<int> sdchR0();\n" <>
          "std::vector<DSlotOpt> optp0();\n" <> "std::vector<std::vector<int>> sdslR0();\n",
          ""
        ] <>
(* chunk helpers of the oversized builders — the bin-packer may put a builder's helpers in a
   different unit than its assembler, so these must be declared here, not per-unit. *)
        chunkDecls;
    (* component-table init: comp[base][mu] = <MPoly builder>, skipping structural zeros. *)
    compInit =
      StringJoin @
        KeyValueMap[
          Function[{base, comps},
            StringJoin @
              MapIndexed[
                Function[{s, mu},
                  If[s === "env.zero()",
                    "",
                    "  comp[" <> str[base] <> "][" <> str[mu[[1]] - 1] <> "] = " <> s <> ";\n"]],
                comps]],
          compCpp];
    main =
      StringJoin[
        Flatten[
          {
            "int main(int argc, char** argv){\n",
            "  std::string decor = \"static inline\"; std::string hns = \"" <> nsInner <> "\";\n",
            "  for(int a=1;a<argc;++a){ std::string s=argv[a]; if(s==\"-d\"&&a+1<argc) decor=argv[++a]; else if(s==\"-n\"&&a+1<argc) hns=argv[++a]; }\n",
            "  const int nsym = " <> str[nsym] <> ";\n",
(* units is emitted BEFORE the LorentzEnv because the env binds both nsym and the unit groups;
   comp/atomDen and the trace entry points are then built through `env`. *)
            "  std::vector<std::vector<int>> units = {" <> StringRiffle[("{" <> StringRiffle[str /@ #, ","] <> "}")& /@ unitG, ","] <> "};\n",
            "  LorentzEnv env(nsym, units);\n",
            "  std::vector<std::array<MPoly,4>> comp(" <> str[maxBase + 1] <> ", {env.zero(),env.zero(),env.zero(),env.zero()});\n",
            compInit,
            "  std::vector<std::string> symNames = {" <> StringRiffle[("\"" <> # <> "\"")& /@ symNames, ","] <> "};\n",
(* the DISTINCT-trace tables (see the global sub-term dedup). Flat, indexed by trace id: a plain
   trace has an empty chain (contract via sdn[k]), a dressed (structural) one uses sdch[k]/sdsl[k] via
   numeric_value_dressed_netval_mp. sdch/sdsl are emitted only when the kernel has dressed nets. *)
            "  std::vector<DiracNet> sdn = sdn0();\n",
            "  std::vector<NetVal> sln = sln0();\n",
            If[hasDressed,
(* rebuild the per-sub-term chain/slot tables from the interned pools (chp/optp) + index arrays
   (sdchR/sdslR) — O(nSub), reproduces the full sdch/sdsl exactly, so the trace lambda is untouched. *)
              "  std::vector<std::vector<DChainTok>> chp = chp0(); std::vector<int> sdchR = sdchR0();\n" <>
              "  std::vector<DSlotOpt> optp = optp0(); std::vector<std::vector<int>> sdslR = sdslR0();\n" <>
              "  const size_t NSD = sdchR.size();\n" <>
              "  std::vector<std::vector<DChainTok>> sdch(NSD); std::vector<std::vector<DSlot>> sdsl(NSD);\n" <>
              "  for(size_t k=0;k<NSD;++k){ sdch[k]=chp[sdchR[k]]; sdsl[k].reserve(sdslR[k].size());\n" <>
              "    for(int oi: sdslR[k]) sdsl[k].push_back(DSlot{optp[oi]}); }\n",
              ""],
(* per net: which traces it references, and with what scalar (sub-terms sharing a trace have
   already been merged, and zero sums dropped, at codegen time) *)
(* sidx/dsc HASH-CONSED. Written out in full these two tables dominate the main TU: on the
   four-quark Fierz gate they were 0.29 MB and 2.14 MB of literals for 3470 nets, and since the
   non-dressed main TU compiles at -O2 (see mainOpt) a single multi-megabyte braced-init costs
   MINUTES — measured >600 s at -O2 vs 6 s at -O0 on the same file, while every net-builder unit
   took 0.65 s. $ntDefChunk already solved this for the net-builder units; these data tables were
   never covered by it.
   The redundancy is extreme because sub-terms sharing a trace share their scalars: 73144 Cx
   literals over only 601 DISTINCT values, and 3470 index rows over 366 distinct. So dedupe rows
   for sidx, and for dsc dedupe BOTH levels (values, then the rows of value-indices — row dedup
   alone leaves ~1 MB, since the rows differ while their entries repeat). The runtime rebuild
   below is O(nets) and reproduces sidx/dsc EXACTLY as before, so fold_nets is untouched. *)
            With[{
              idxRows =
                Function[nt,
                    nt[[All, 1]]
                  ] /@ netTerms,
              (* LEVER (b): the per-sub-term dressing monomials (position 2), row-deduped like sidx *)
              drRows =
                Function[nt,
                    nt[[All, 2]]
                  ] /@ netTerms,
              scaRows =
                Function[nt,
                    nt[[All, 3]]
                  ] /@ netTerms},
              With[{distinctIdxRows = DeleteDuplicates[idxRows], distinctScalars = DeleteDuplicates[Flatten[scaRows]], distinctDressMonos = DeleteDuplicates[Flatten[drRows, 1]]},
                With[{idxPos = AssociationThread[distinctIdxRows -> Range[Length[distinctIdxRows]] - 1], valPos = AssociationThread[distinctScalars -> Range[Length[distinctScalars]] - 1], drValPos = AssociationThread[distinctDressMonos -> Range[Length[distinctDressMonos]] - 1]},
                  With[{scaIdxRows = Map[valPos, scaRows, {2}], drIdxRows = Map[drValPos, drRows, {2}]},
                    With[{distinctScalarRows = DeleteDuplicates[scaIdxRows], distinctDressRows = DeleteDuplicates[drIdxRows]},
                      With[{scaPos = AssociationThread[distinctScalarRows -> Range[Length[distinctScalarRows]] - 1], drRowPos = AssociationThread[distinctDressRows -> Range[Length[distinctDressRows]] - 1]},
                        StringJoin[
                          "  const size_t NNET = " <> str[Length[netTerms]] <> ";\n",
                          "  std::vector<std::vector<int>> sidxU = {" <> StringRiffle[("{" <> StringRiffle[str /@ #, ","] <> "}")& /@ distinctIdxRows, ","] <> "};\n",
                          "  std::vector<int> sidxR = {" <> StringRiffle[str /@ (idxPos /@ idxRows), ","] <> "};\n",
                          "  std::vector<std::vector<int>> sidx(NNET);\n",
                          "  for(size_t i=0;i<NNET;++i) sidx[i]=sidxU[sidxR[i]];\n",
                          "  std::vector<Cx> dscV = {" <>
                            StringRiffle[
                              Function[s,
                                  "Cx{" <> cppNum[Re[s]] <> "," <> cppNum[Im[s]] <> "}"
                                ] /@ distinctScalars,
                              ","
                            ] <> "};\n",
                          "  std::vector<std::vector<int>> dscU = {" <> StringRiffle[("{" <> StringRiffle[str /@ #, ","] <> "}")& /@ distinctScalarRows, ","] <> "};\n",
                          "  std::vector<int> dscR = {" <> StringRiffle[str /@ (scaPos /@ scaIdxRows), ","] <> "};\n",
                          "  std::vector<std::vector<Cx>> dsc(NNET);\n",
                          "  for(size_t i=0;i<NNET;++i){ const auto& r=dscU[dscR[i]]; dsc[i].reserve(r.size());\n",
                          "    for(int k: r) dsc[i].push_back(dscV[k]); }\n",
(* LEVER (b): the per-sub-term dressing monomials sdr[i][j], deduped on BOTH levels exactly like dsc — a
   dense dressed flow has thousands of sub-terms but only a handful of DISTINCT dressing monomials (sdrV)
   and few distinct index-rows (sdrU), so a flat braced-init would be huge (measured: 115 KB single row on
   za3_147, +37 s compile) while this stays a few KB. The dressed phase-B fold routes each sub-term's scaled
   MPoly trace into its DPoly channel sdr[i][j] (empty monomial = undressed). Emitted only for dressed flows. *)
                          If[hasDressed,
                            "  std::vector<DMono> sdrV = {" <>
                              StringRiffle[("{" <> StringRiffle[str /@ #, ","] <> "}")& /@ distinctDressMonos, ","] <> "};\n" <>
                            "  std::vector<std::vector<int>> sdrU = {" <>
                              StringRiffle[("{" <> StringRiffle[str /@ #, ","] <> "}")& /@ distinctDressRows, ","] <> "};\n" <>
                            "  std::vector<int> sdrR = {" <> StringRiffle[str /@ (drRowPos /@ drIdxRows), ","] <> "};\n" <>
                            "  std::vector<std::vector<DMono>> sdr(NNET);\n" <>
                            "  for(size_t i=0;i<NNET;++i){ const auto& r=sdrU[sdrR[i]]; sdr[i].reserve(r.size());\n" <>
                            "    for(int k: r) sdr[i].push_back(sdrV[k]); }\n",
                            ""]]]]]]]],
(* the projector atom denominators are keyed by ATOM ID (e.inv/e.invS), not by position, and each
   id is filled idempotently — and the distinct traces cover every lnet that occurs. So scanning
   the deduped table gives the same atomDen as scanning every occurrence did. *)
            "  auto atomDen = env.collect_atom_denoms(sln, comp);\n",
            "  for(auto &a: atomDen) a = reduce_units(a, units);  // bare-loop k^2 -> monomial l1^2 -> cancels\n",
            "  const bool ntprof = (std::getenv(\"NT_GEN_PROFILE\")!=nullptr);\n",
            "  unsigned hw=std::thread::hardware_concurrency(); if(!hw)hw=4u;\n",
            "  if(const char* mw=std::getenv(\"NT_GEN_MAXW\")){int v=std::atoi(mw); if(v>0&&(unsigned)v<hw)hw=(unsigned)v;}\n",
(* SEPARATE worker count for phase B. The two phases have very different memory profiles per
   worker, so one knob cannot tune both. Phase A's transient is `hw` concurrent contractions, but
   it is trimmed away afterwards. Phase B's is `hw` concurrent RECOMPUTES of the uncached traces —
   and those are the SINGLETONS, which are the heavy ones (Codegen.m orders by descending refcount,
   and a trace that recurs is a simple structure while a unique one is complex: on ZAAqbq1 the 2550
   reused traces are 20.1 MB total, ~8 KB each, while the 1164 singletons could not be cached at
   all inside 10 GB). Those recomputes land on top of the live window, so phase B can need FEWER
   workers than phase A even though it is the cheaper phase in CPU terms. Defaults to `hw`. *)
            "  unsigned hwB=hw; if(const char* mb=std::getenv(\"NT_GEN_MAXW_B\")){int v=std::atoi(mb); if(v>0)hwB=(unsigned)v;}\n",
            With[{
              (* LEVER (b): the trace table is PLAIN MPoly for BOTH paths now — a dressed sub-term's
                        structural trace is a plain MPoly (numeric_value_dressed_netval_mp) and its dressing
                        rides the per-sub-term scalar (dsc) + monomial (sdr), assembled into a DPoly only in
                        the phase-B fold. So phase A caches MPoly either way. *)
              PT = "MPoly"},
              StringJoin[
                {
                  "  const long NSUB = " <> str[nSub] <> ";\n",
(* how many traces are RESIDENT. Default: the reused ones (refCount >= 2), which the codegen-time
   ordering puts first — a singleton is contracted once whether cached or not, so caching it is
   pure RAM for no saving. NT_GEN_MEMO_MAX overrides either way (clamped to [0, NSUB]): lower it
   when memory is tight (the RAM lever — the dense flows are memory-bound before they are
   compute-bound), raise it to NSUB to put the singletons in phase A too, which costs their RAM
   but gives phase A the whole work list to balance.
   DRESSED (collected-slot) flows default to NSUB: even after lever (b) dedups the dressing variants of a
   concrete trace (so nReused is no longer ~0), the plain-MPoly traces are individually SMALL and phase B
   is parallel over NETS not traces — leaving the singletons to phase B lets the one dominant net serialise
   thousands of contractions. Caching them all (nSub) is RAM-cheap now the trace table has no dressing
   dimension, and lets phase A contract them over its flat W-parallel work list. Memory-bound dressed flows
   (ZAAqbq) dial it back with NT_GEN_MEMO_MAX. *)
                  "  long nCache = " <> str[If[hasDressed, nSub, nReused]] <> ";\n",
                  "  if(const char* mm=std::getenv(\"NT_GEN_MEMO_MAX\")){ long v=std::atol(mm); if(v>=0) nCache=std::min<long>(v,NSUB); }\n",
                  "  auto trace=[&](int k)->" <> PT <> "{\n",
                  If[hasDressed,
                    (* structural trace → plain MPoly: a non-slot sub-term contracts via sdn[k]/sln[k],
                                 a collected one via the dressing-free _mp variant (dressing stripped at codegen). *)
                    "    return sdch[k].empty()\n" <> "      ? env.numeric_value_netval(sdn[k], sln[k], comp, atomDen)\n" <> "      : env.numeric_value_dressed_netval_mp(sdch[k], sdsl[k], sln[k], comp, atomDen);\n",
                    "    return env.numeric_value_netval(sdn[k], sln[k], comp, atomDen);\n"],
                  "  };\n",
                  (* PHASE A — contract each distinct trace once, parallel over a FLAT work list (numeric/trace_fold.hpp). *)
                  "  auto tA=std::chrono::steady_clock::now();\n",
                  "  std::vector<" <> PT <> "> T = env.contract_traces<" <> PT <> ">(nCache, hw, trace);\n",
                  "  if(ntprof){ std::size_t tb=0; for(auto &p: T) tb+=poly_bytes(p);\n",
                  "    std::fprintf(stderr,\"[num] phase A: %ld distinct traces, %ld cached, table %.1f MB, %.1f s (W=%u)\\n\",\n",
                  "      NSUB, nCache, tb/1048576.0, std::chrono::duration<double>(std::chrono::steady_clock::now()-tA).count(), hw); }\n",
(* RELEASE PHASE A'S ARENA. Phase A contracts `hw` traces CONCURRENTLY, and a single dense 4-point
   contraction transiently allocates ~1 GB — so its working set is ~hw GB even though the table it
   leaves behind is 20 MB. glibc frees that into the per-thread arenas but does not munmap it, so
   RSS stays at the phase-A high-water mark and phase B starts from a multi-GB floor instead of
   from the table. Measured on ZAAqbq1 at W=6: RSS ~7.4 GB after phase A against a 20.1 MB table,
   which is what then pushed phase B over a 10 GB cap. malloc_trim gives it back to the OS.
   Costs milliseconds, once. *)
                  "#if defined(__GLIBC__)\n",
                  "  { const double rssPre = ntRssMB(); malloc_trim(0);\n",
                  "    if(ntprof) std::fprintf(stderr,\"[num] arena trim after phase A: RSS %.0f -> %.0f MB\\n\", rssPre, ntRssMB()); }\n",
                  "#endif\n",
(* PHASE B is NOT emitted here — it is fused into the group/lowering loop below (search
   fold_groups_streaming). It used to be `mp = env.fold_nets(...)`, one fully-expanded polynomial per
   net, ALL of them returned and then held for the rest of main() while the group loop summed them
   into a second full set of per-group accumulators. On the dense 4-point flows that is 20+ GB (488
   nets x ~41 MB) against a 20 MB trace table — the generator was OOM-killed before it could emit.
   Nothing is revisited (each net polynomial is written once, read once by its group, then dead), so
   the fix is to consume it as a stream: fold each group's nets on demand, lower the group, free it.
   Phase B therefore needs `groups`, `colv`, `g` and `realOnly`, which are only declared further
   down — hence the move. `tB` still starts here so the reported phase-B time is comparable. *)
                  "  auto tB=std::chrono::steady_clock::now();\n"}]],
(* one SUNEnv per distinct group rank appearing in the colour nets (colFacG emits `sun<n>.` factors),
   so the rank is written once. Scan the assembled net strings for the `sun<digits>.` tokens. *)
            StringJoin["  SUNEnv sun" <> # <> "(" <> # <> ");\n"& /@ DeleteDuplicates @ Flatten @ StringCases[colourNets, "sun" ~~ r : DigitCharacter.. ~~ "." :> r]],
(* colnets HASH-CONSED for the same reason (0.57 MB -> 0.115 MB on the Fierz gate: 3470 rows,
   719 distinct). Two nets that differ only in their Lorentz/Dirac part share a colour net, so
   the duplication is structural. This also collapses the sun_value_cx calls to the distinct
   nets — 4.8x fewer on that flow — which is a RUN saving on top of the compile one. *)
            With[{uCol = DeleteDuplicates[colourNets]},
              With[{colPos = AssociationThread[uCol -> Range[Length[uCol]] - 1]},
                StringJoin["  std::vector<SUNNet> colnetsU = {" <> StringRiffle[uCol, ","] <> "};\n", "  std::vector<int> colR = {" <> StringRiffle[str /@ (colPos /@ colourNets), ","] <> "};\n", "  std::vector<numtracer::Cx> colvU(colnetsU.size());\n", "  for(size_t i=0;i<colnetsU.size();++i) colvU[i]=sun_value_cx(colnetsU[i]);\n", "  std::vector<numtracer::Cx> colv(" <> str[nNet] <> ");\n", "  for(int i=0;i<" <> str[nNet] <> ";++i) colv[i]=colvU[colR[i]];\n"]
              ]],
            "  std::vector<std::vector<int>> groups = {" <> StringRiffle[("{" <> StringRiffle[str /@ #, ","] <> "}")& /@ groups, ","] <> "};\n",
            "  GlobalEnv g;\n",
            "  std::vector<GenProg> progs;\n",
(* realOnly[gi]: this group's dressing coeff is real, so only Re(trace) is consumed -> emit a
   double trace (no dead imaginary half). Defaults to all-0 (full complex) if not supplied. *)
            "  std::vector<int> realOnly = {" <>
              StringRiffle[
                If[Length[realOnlyG] === nGrp,
                  (
                      If[TrueQ[#],
                        "1",
                        "0"]
                    )& /@ realOnlyG,
                  Table["0", {nGrp}]],
                ","
              ] <> "};\n",
(* CrossTraceCSE: accumulate every group's polynomial first, then lower them all through ONE
   shared CSE builder (to_genprog_fused) instead of one independent program per trace. Measured on
   ZAqbq1_147 Mq-in: 30,547 shared SSA instrs vs 47,558 independent (0.64x), lowering cost
   unchanged. See network::FusedProg. *)
(* PHASE B + group accumulation + lowering, FUSED into one streaming pass (numeric/trace_fold.hpp's
   fold_groups_streaming). `groups` partitions the nets, so a net's folded polynomial is needed by
   exactly one group and can die the moment that group has absorbed it. At most `gwin` group
   accumulators plus `hw` in-flight net polynomials are ever live, against every net polynomial AND
   every group accumulator before. NT_GEN_GROUP_WINDOW is the dial; gwin == nGrp reproduces the old
   residency exactly, so the pre-streaming behaviour stays reachable for A/B.

   Bit-identical to what it replaces: each group still left-folds its members in group order over the
   same fold_net results, and the sink runs on the calling thread for gi = 0,1,2,... ascending — which
   is what keeps GlobalEnv intern order and the shared CSE instruction stream unchanged. The scale is
   passed per branch rather than shared so each keeps its exact expression (poly*constant here vs
   scale_trace's constant*poly). *)
            "  long gwin = numtracer::numeric::net_window((long)sidx.size(), hwB);\n",
            "  if(ntprof) numtracer::numeric::check_group_partition(groups, " <> str[nNet] <> ");\n",
(* LEVER (b): the dressed fold reads a PLAIN-MPoly trace table T + the per-sub-term dressing monomials sdr,
   building the per-net DPoly channel-by-channel (fold_groups_streaming_dressed). The colour scale (scaleCx)
   and the sink are unchanged — only the per-net fold's trace type differs, so the group DPoly the sink
   receives is value-identical to the pre-lever-(b) DPoly-trace-table path. *)
            Which[
              crossCSE && hasDressed,
                "  FusedStream fstream(g, realOnly);\n" <> "  env.fold_groups_streaming_dressed(sidx, dsc, sdr, groups, T, nCache, hwB, gwin, trace,\n" <> "    [&](int d, DPoly &&m){ return scaleCx(m, colv[d]); },\n" <> "    [&](size_t, DPoly &&acc){ fstream.add(acc); });\n" <> "  std::vector<FusedProg> fused = fstream.finish();\n",
              crossCSE,
                "  FusedStream fstream(g, realOnly);\n" <> "  env.fold_groups_streaming<MPoly>(sidx, dsc, groups, T, nCache, hwB, gwin, trace,\n" <> "    [&](int d, MPoly &&m){ return m*env.constant(colv[d]); },\n" <> "    [&](size_t, MPoly &&acc){ fstream.add(acc); });\n" <> "  std::vector<FusedProg> fused = fstream.finish();\n",
              hasDressed,
                "  env.fold_groups_streaming_dressed(sidx, dsc, sdr, groups, T, nCache, hwB, gwin, trace,\n" <> "    [&](int d, DPoly &&m){ return scaleCx(m, colv[d]); },\n" <> "    [&](size_t gi, DPoly &&acc){ progs.push_back(to_genprog(acc, g, realOnly[gi]!=0)); });\n",
              True,
                "  env.fold_groups_streaming<MPoly>(sidx, dsc, groups, T, nCache, hwB, gwin, trace,\n" <> "    [&](int d, MPoly &&m){ return m*env.constant(colv[d]); },\n" <> "    [&](size_t gi, MPoly &&acc){ progs.push_back(to_genprog(acc, g, realOnly[gi]!=0)); });\n"
            ],
            "  if(ntprof) std::fprintf(stderr,\"[num] phase B+lower: %d nets in %d groups, window %ld, %.1f s (W=%u)\\n\", " <> str[nNet] <> ", " <> str[nGrp] <> ", gwin, std::chrono::duration<double>(std::chrono::steady_clock::now()-tB).count(), hwB);\n",
(* the trace table is dead once every group has folded; emission below only needs the lowered
   instruction streams, which are orders of magnitude smaller. *)
            (* LEVER (b): the trace table T is plain MPoly for BOTH paths now. *)
            "  { std::vector<MPoly> dead; T.swap(dead); }\n",
            "  FillFormulas fm;\n",
            "  fm.var = [](int id)->std::string{\n",
            Table["    if(id==" <> str[i - 1] <> ") return \"" <> varFill[[i]] <> "\";\n", {i, 1, Length[varFill]}],
            "    return \"0.0\"; };\n",
            "  fm.inv = [&](int id)->std::string{ return \"1.0/(\" + mpoly_to_cpp(atomDen[(size_t)id], symNames) + \")\"; };\n",
(* dressing fill (kind-2 `dress` leaves): the kernel BODY evaluates each dressing atom (where the
   regulators REG::* and the interpolator parameters are in scope) into `dr_<id>` and passes the
   VALUE to fill(); fill just stores it. So fm.dress returns the passed-in argument name. *)
            If[hasDressed,
              "  fm.dress = [](int id)->std::string{ return \"dr_\" + std::to_string(id); };\n",
              ""],
            "  std::cout << \"// GENERATED by gen_" <> nsInner <> ".cpp — do not edit.\\n\";\n",
            "  std::cout << \"#pragma once\\n#include <cmath>\\n" <>
              If[complexQ,
                "#include <complex>\\n",
                ""
              ] <> "namespace " <> kns <> " { namespace \" << hns << \" {\\n\";\n",
            "  std::cout << \"template<int N> \" << decor << \" double powr(double x){ double r=1.0; for(int i=0;i<N;++i) r*=x; return r; }\\n\";\n",
            "  emit_env_layout(std::cout, g);\n",
            "  std::cout << \"static inline constexpr int nenv = \" << g.syms.size() << \";\\n\";\n",
            "  emit_fill(std::cout, g, \"fill\", \"" <> fillArgSig <> "\", fm, decor);\n",
(* TRACE-BODY DEDUP: the grouping key is the dressing COEFFICIENT (diagData), finer than the trace
   STRUCTURE — so flows with many Feynman graphs that share a kinematic trace but differ only in
   their dressing/coupling coefficient (e.g. the σ-Yukawa hSigL: 503 groups, only 157 distinct trace
   bodies) emit hundreds of byte-identical trN. Render each trace, key on the name-independent body,
   and emit a duplicate as a one-line forwarder `trN(f){ return trK(f); }` (canonical K = first with
   that body). The body is computed identically once per call site regardless (GCC CSEs the inlined
   identical traces — see ZA4 fusion notes), so this is a pure SOURCE-SIZE/compile win, runtime-neutral.
   Flows with no shared trace structure (ZAqbq1/4/7_147: 108/108 distinct) never hit `seen` ⇒ the
   emitted bytes are unchanged. *)
(* fused: ONE trace_all(f, t[]) instead of nGrp trN(). The body-dedup below is meaningless then
   (there is a single body), and the kernel reads tarr[i] rather than calling tr_i. *)
            If[crossCSE,
              "  emit_cpp_fused(std::cout, fused, \"trace_all\", decor);\n",
              "  { std::unordered_map<std::string,std::string> seen; seen.reserve((size_t)" <> str[nGrp] <> ");\n" <> "    for(int i=0;i<" <> str[nGrp] <> ";++i){\n" <> "      std::ostringstream os; emit_cpp(os, progs[i], \"tr\"+std::to_string(i), decor);\n" <> "      std::string s = os.str(); std::string body = s.substr(s.find('{'));\n" <> "      auto it = seen.find(body);\n" <> "      if(it==seen.end()){ seen.emplace(std::move(body), \"tr\"+std::to_string(i)); std::cout << s; }\n" <> "      else { const char* rt = (s.find(\"std::complex<double> tr\")!=std::string::npos) ? \" std::complex<double> \" : \" double \";\n" <> "        std::cout << decor << rt << \"tr\" << i << \"(const double *f) { return \" << it->second << \"(f); }\\n\"; } } }\n"
            ],
            "  std::cout << \"}} // namespace " <> kns <> "::\" << hns << \"\\n\";\n  return 0;\n}\n"}]];
    {pre, units, decl, main}];

(* ---- whole kernel (boilerplate delegated to FunKit) ------------------------- *)
(* regulator wrappers, prefixed with the user-chosen decorator (default plain "static inline";
   pass e.g. "Decorator" -> "static __host__ __device__ inline" for CUDA-callable kernels).
   Emitted ONLY under "RegulatorTemplate" -> True (the fRG/DiFfRG shape); see ntKernelClass. *)

privDefs[decor_] :=
  StringRiffle[(decor <> " auto " <> # <> "(const auto &k2, const auto &p2) { return REG::" <> # <> "(k2, p2); }")& /@ {"RB", "RF", "RBdot", "RFdot", "dq2RB", "dq2RF"}, "\n"];

(* per-trace real/imag accessors for the "RePart" re/im split: a kernel whose VALUE is real but which
   has a COMPLEX trace is assembled as Σ[Re(c)·tr.real() − Im(c)·tr.imag()] — PURE double arithmetic, so
   no complex type survives into the kernel (avoids mixing the trace's std::complex with the support
   complex, and keeps device code real). Overloads pass a real trace straight through (im → 0). *)

ntReImDefs[decor_] :=
  StringRiffle[{decor <> " double ntRe(double x) { return x; }", "template <class T> " <> decor <> " double ntRe(const T &z) { return z.real(); }", decor <> " double ntIm(double) { return 0.0; }", "template <class T> " <> decor <> " double ntIm(const T &z) { return z.imag(); }"}, "\n"];

(* ---- the kernel CLASS. -----------------------------------------------------------------------
   Two shapes, one code path (it is emitted twice — once up front, once after the real probe
   re-lowers the integrand — and the two must not drift):

     regTemplate = False (DEFAULT, the general emission): a PLAIN class. NumTracer emits the kernel
       and nothing else; a flow whose dressing rules mention RB/RF/RBdot/RFdot/dq2RB/dq2RF emits
       UNQUALIFIED calls to them, which the consumer supplies (e.g. via "ExtraIncludes" -> {hdr},
       or from the "RuntimeInclude" header). Unqualified lookup runs class -> kernel namespace ->
       global, so free functions at global scope are found from any "KernelNamespace".

     regTemplate = True (the fRG/DiFfRG shape): `template<typename REG> class ...` plus the six
       private wrappers forwarding to REG::. DiFfRG's scaffold forward-declares the kernel as a
       template and instantiates it as KERNEL<Regulator>, so MakeNTKernelDiFfRG needs this.

   regAlias adds `using Regulator = REG;` (DiFfRG reads it back off the kernel class) and only makes
   sense with the template, which is why it implies it at the call sites below. *)

ntKernelClass[name_, members_List, decor_, regTemplate_, regAlias_, extraPriv_List] :=
  FunKit`MakeCppClass[
    Sequence @@
      If[TrueQ[regTemplate],
        {"TemplateTypes" -> {"REG"}},
        {}],
    "Name" -> name,
    "MembersPublic" ->
      If[TrueQ[regAlias],
        Prepend[members, "using Regulator = REG;"],
        members],
    "MembersPrivate" ->
      If[TrueQ[regTemplate],
        Join[{privDefs[decor]}, extraPriv],
        extraPriv]];

(* ---- emit helpers: support `using`s, namespace wrapping, runtime include -------------------
   A generated kernel pulls its math helpers (complex, powr, pow, sqrt, fma) from a configurable
   SUPPORT namespace and is wrapped in a configurable HOST namespace. The defaults (numtracer /
   numtracer/codegen/runtime.hpp) make the emitted code self-contained against NumTracer's own
   headers; a consumer that already provides equivalents (e.g. DiFfRG) points the codegen at
   them via "SupportNamespace" / "KernelNamespace" / "RuntimeInclude". *)

ntSupportUsings[sns_] :=
  DeleteDuplicates[{"using namespace " <> sns <> ";", "using namespace " <> sns <> "::compute;", "using namespace numtracer;"}];

(* the loop-independent `constant(p, k, dressings...)` function. DiFfRG flat-adds its return to the
   integral (constExpr second arg of MakeKernel); NumTracer emits the same. A 0 constant needs no
   namespace usings (and stays byte-identical to the pre-Constant emission); a nonzero expr may call
   compute helpers (powr/pow/…) or the support API, so it gets the same usings the kernel body does. *)

ntConstFn[constExpr_, decor_, constParams_, sns_] := FunKit`MakeCppFunction[
    constExpr,
    "Name" -> "constant",
    "Prefix" -> decor,
    "Return" -> "auto",
    "CodeParser" -> "Cpp",
    "Parameters" -> constParams,
    "Body" ->
      If[MatchQ[constExpr, 0 | 0.],
        "",
        StringRiffle[ntSupportUsings[sns], "\n"]]];

ntWrapBody[kns_, classStr_, name_] := If[kns === None || kns === "",
    {classStr},
    {"namespace " <> kns <> "\n{", classStr, "}", "using " <> kns <> "::" <> name <> ";"}];

ntRuntimeIncludes[runInc_] := If[runInc === None || runInc === "",
    {},
    {runInc}];

(* the dressing-parameter type: Automatic -> `const auto&` (fully generic, self-contained);
   else the given concrete type string (e.g. a consumer's interpolator type). *)

ntDressType[dressTy_] := If[dressTy === Automatic,
    "auto",
    dressTy];

(* ---- numeric matrix-product kernel, the MakeNTKernel "Numeric" path: emit the build-time
        generator program + the kernel header that calls its output. ---------------------- *)

Options[mkGenerateKernel] =
  {
    "Name" -> "nt_inv_kernel",
    "Namespace" -> Automatic,
    "Dressings" -> {},
    "ScalarParams" -> {},
    "ADParams" -> {},
    "IncludeDir" -> Automatic,
    "RunGenerator" -> True,
    "FullParallel" -> False,
    "AngleDefs" -> {},
    "CrossTraceCSE" -> False,
    "GlobalCollect" -> True,
    (* default ON: groups diagrams by dressing coeff + folds colour
numerically so COEN collects the factored coefficient across diagrams (FORM-style); ~30% runtime
cut on the dense quark-loop vertices, generation cost unchanged. Pass False to opt out per flow. *)
    "NumericContract" -> False,
    "Components" -> Automatic,
    "SymbolDefs" -> <||>,
    "Decorator" -> "static inline",
    "RuntimeInclude" -> "numtracer/codegen/runtime.hpp",
    "ExtraIncludes" -> {},
    "KernelNamespace" -> "numtracer_kernels",
    "SupportNamespace" -> "numtracer",
    "DressingType" -> Automatic,
    "RegulatorTemplate" -> False,
    "RegulatorAlias" -> False,
    "RealProbe" -> True,
    "PruneRealTraces" -> False,
    "Constant" -> 0.};

(* The dressing-collection path is driven by the `ntDressedNum` tokens NumTrace emits under
   "DressingCollection" -> True (no MakeNTKernel option needed): mkGenerateKernel auto-detects them and
   routes to the DPoly generator branch. *)
(* "Constant" -> expr (default 0.): the loop-INDEPENDENT piece that DiFfRG flat-adds to the
   integral — emitted verbatim as the body of `constant(p, k, dressings...)`, exactly like the
   constExpr second argument of DiFfRG's MakeKernel. A plain Mathematica expression in p/k and the
   dressing names (e.g. ZA[p] -> ZA(p)); NOT an NTKernel. Left at 0. the constant returns 0.. *)
(* "PruneRealTraces" -> False (default OFF): emit a `double` trace for any diagram group whose dressing
   coefficient is real (only Re(trace) is consumed), skipping the dead imaginary half. CORRECT for the
   kernel in every verdict, but HAZARDOUS with the probe: the probe runs on the generated traces and
   verifies Im(integrand)≈0; dropping a real-coeff group's imaginary half removes a term the probe
   needs, which can leave a non-cancelling residual and misclassify a real flow as Complex (blocking
   the double-kernel emission). Safe to enable only when the probe is off OR after per-flow validation,
   or once a post-probe trace regeneration is wired. Left off pending measurement (GCC often already
   DCEs the dead half of small inlined traces). *)
(* "RealProbe" -> True (default): when the syntactic complexQ trips (some diagram coeff carries an `i`),
   compile+run a probe over the JUST-generated real traces to test whether Im(integrand) actually
   vanishes (projector-i × colour-i often cancel to a real value that Mathematica can't see through the
   opaque trace symbols). If it vanishes, re-emit a REAL (double) kernel — losslessly. Set False to skip
   the probe and always keep the complex+consumer-Re path. *)
(* "NumericContract" -> True: the numeric (matrix-product) backend (task #22) — contract each diagram
   numerically via et/numeric (γ/metric/projector numeric, momenta symbolic), no sp-invariant
   reduce/rebase/ibp. "Components" -> <|mom -> {e0,e1,e2,e3}, ...|> gives each momentum's 4 components
   as expressions (partially numeric / partially symbolic); Automatic falls back to the kinematic
   frame polynomialised. "SymbolDefs" -> <|sym -> expr|> gives the C++ fill for any DERIVED symbol
   (e.g. sin1 -> Sqrt[1-cos1^2]); plain free symbols are taken to be kernel arguments. *)
(* Standalone-output options (defaults make the emitted code self-contained against NumTracer's
   own headers, with no mention of any downstream consumer). A consumer that supplies its own
   support API points the codegen at it via these:
     "RuntimeInclude" -> "<hdr>" | None : the support header #included first, providing `complex`
        and `compute::{powr,pow,sqrt,fma}` (default numtracer/codegen/runtime.hpp; None to omit).
     "ExtraIncludes" -> {"a.hpp", ...}  : extra #includes prepended ahead of everything.
     "KernelNamespace" -> "ns" | None   : namespace wrapping the kernel class AND the generated
        trace functions (default "numtracer_kernels"; None emits at the includer's scope).
     "SupportNamespace" -> "ns"         : where `complex`/`compute` are looked up via `using`
        (default "numtracer").
     "DressingType" -> Automatic | "T"  : dressing-parameter type; Automatic emits `const auto&`
        (fully generic), or give a concrete type string.
     "RegulatorTemplate" -> False       : by default the kernel is a PLAIN class. A flow whose
        dressing rules use the regulators emits unqualified RB/RF/RBdot/RFdot/dq2RB/dq2RF calls,
        which the CONSUMER supplies — put them at global (or kernel-namespace) scope in a header and
        pull it in with "ExtraIncludes" -> {"my_regulators.hpp"} (or from the "RuntimeInclude"
        header). True restores the fRG/DiFfRG shape: `template<typename REG> class ...` plus private
        wrappers forwarding to REG::. Implied by "RegulatorAlias".
     "RegulatorAlias" -> False          : emit `using Regulator = REG;` in the class (DiFfRG reads
        the regulator type back off the kernel). Implies "RegulatorTemplate" -> True. *)
(* "RunGenerator" -> False: emit the generator sources only, without compiling/running them
   (the committed traces header is left untouched). *)
(* "Decorator" -> "<prefix>": the function prefix on EVERY emitted function — kernel/constant,
   the regulator wrappers, and (via the generator's runtime `-d` flag) fill/trN/powr in the
   straight-line header — e.g. "static __host__ __device__ inline" makes the whole kernel
   CUDA-device-callable. Default keeps the emitted bytes identical (kernel md5 invariant). *)
(* "CrossTraceCSE" -> True: lower all diagram trace polynomials through ONE shared CSE program
   (emitted as trace_all(f, t[]); the kernel fills t[] once and reads t[d]) so subexpressions are
   shared ACROSS traces — a fused kernel. Pays off when many traces share intermediates; default
   False keeps one independent trN() per trace.
   Works on COMPLEX flows: t[] is typed by the emitted `trace_all_t` (std::complex<double> iff some
   trace lowered complex), so nothing is truncated.
   Measured ZAqbq1_147 Mq-in (108 traces, 54 complex): 30,547 shared SSA instrs vs 47,558
   independent = 0.64x, lowering cost unchanged. The cross-diagram monomial duplication is 3.54x but
   that is NOT the attainable factor — each monomial occurrence still needs its own accumulate into
   its own trace; only the products are shared. Watch IPC as well as instruction count: this collapses
   N small functions into one very large basic block, which can spill (ZA "Route-B" shrank code 2.2x
   and regressed runtime 31%). *)
(* "AngleDefs" -> {sym -> expr, ...}: kinematic angle symbols the dressing keeps SYMBOLIC
   (e.g. cosl1p2 -> (-cos1 + Sqrt[3-3 cos1^2] cos2)/2). Emitted ONCE as `const double sym = ...;`
   in the kernel body so a shared sub-expression (the sqrt) is computed once rather than inlined
   per occurrence. *)
(* "FullParallel" -> True passes `-p` to the generator: heavy nets are reduced/rebased CONCURRENTLY
   (faster codegen, higher peak RAM) instead of one-at-a-time. The emitted kernel is identical. *)

mkGenerateKernel::genfail = "Generator compile/run failed: `1`";

(* (mkGenerateKernel::crosscseComplex was removed 2026-07-19: CrossTraceCSE now types tarr[] from the
   emitted `trace_all_t`, so a complex flow is no longer truncated and needs no guard.) *)

mkGenerateKernel::emptynets = "Flow `1` produced no generator nets (nets=`2`, groups=`3`) — nothing to emit. Aborting instead of writing a placeholder kernel. (A diagram-build guard such as disconnectmix may have dropped everything, or NumTrace returned no usable diagrams.)";

mkGenerateKernel::scalarleak = "Diagram `1`: a non-numeric factor `2` reached the generator scalar coefficient (a Lorentz tensor that was not resolved by the net builder, e.g. an un-anchored metric contraction). It would be emitted as undeclared C++. Aborting; fix the net build (compileTInv) so the contraction folds numerically.";

(* Locate the NumTracer C++ headers for the generator compile. In-tree the package sits in
   <repo>/numtracer/mathematica/ with the headers in ../include; installed (e.g. under
   $UserBaseDirectory/Applications/NumTracer) the header location is recorded at configure time
   in the CMake-generated sibling NumTracerPaths.m. Overridable per call ("IncludeDir" option)
   or globally (NUMTRACER_INCLUDE_DIR environment variable). *)

resolveIncludeDir::nodir = "Cannot locate the NumTracer C++ headers (tried the NUMTRACER_INCLUDE_DIR environment variable, the in-tree ../include, the installed NumTracerPaths.m record, and ~/.local/share/NumTracer/include). Pass \"IncludeDir\" -> dir to MakeNTKernel.";

resolveIncludeDir[] := Module[{envDir, dir, pathsFile},
    envDir = Environment["NUMTRACER_INCLUDE_DIR"];
    If[StringQ[envDir] && DirectoryQ[envDir],
      Return[envDir]];
    dir = FileNameJoin[{DirectoryName[$NumTracerDirectory], "include"}];
    If[DirectoryQ[dir],
      Return[dir]];
    pathsFile = FileNameJoin[{$NumTracerDirectory, "NumTracerPaths.m"}];
    If[FileExistsQ[pathsFile],
      Get[pathsFile];
      If[StringQ[$NumTracerInstalledIncludeDir] && DirectoryQ[$NumTracerInstalledIncludeDir],
        Return[$NumTracerInstalledIncludeDir]]];
    dir = FileNameJoin[{$HomeDirectory, ".local", "share", "NumTracer", "include"}];
    If[DirectoryQ[dir],
      Return[dir]];
    Message[resolveIncludeDir::nodir];
    Abort[]];

(* Locate the compiled engine library `libNumTracer.a` for the generator LINK. NumTracer ships as a
   compiled static library by default: the heavy engine bodies (numeric contraction / SU(N) fold /
   lowering) are compiled ONCE into this archive, so the emitted generator only parses the declaration
   headers and links the archive instead of re-instantiating and -O2-optimising the whole engine on
   every generation run. Searched: the NT_GEN_LIB environment variable (a full path), the installed
   `<prefix>/lib` sibling of the include dir, and the in-tree `<repo>/numtracer/build` default build
   dir. Returns the path, or $Failed — in which case the generator falls back to a self-contained
   header-only compile (`-DNUMTRACER_HEADER_ONLY=1`), which still works but pays the old compile floor. *)
(* The generator compiles gen_*.cpp against the CURRENT headers in incDir but LINKS the prebuilt
   archive. Those two must come from the same source vintage: the engine's types (MPoly's storage and
   its construction path, SUNNet, the fold buffers) are defined in the headers and compiled into the
   archive, so a header edit that is not followed by a library rebuild is a silent ODR/ABI mismatch —
   the generator RUNS, emits a plausible kernels.hh, and the numbers are wrong. Nothing downstream can
   catch it: kernel.hh only NAMES the traces (s1..sN) and is emitted by the Mathematica layer, so it
   stays byte-identical while every trace body in kernels.hh changes underneath it — which is exactly
   how a stale archive silently corrupted a flow's ZA4 and survived a "kernel.hh unchanged" check.
   mtime is a coarse proxy for provenance, but the failure it guards is silent and total, so err loud:
   a library older than any header it was built from is never trustworthy. *)

MakeNTKernel::stalelib = "The prebuilt engine archive\n  `1`\nis OLDER than the header\n  `2`\nthat the generator will compile against. Linking them mixes two source vintages of the engine's types (MPoly/SUNNet/fold buffers) — an ODR/ABI mismatch that SILENTLY produces wrong traces in kernels.hh while leaving kernel.hh byte-identical. Rebuild the library first, e.g.\n  cmake --build <repo>/numtracer/build --target NumTracer\nthen regenerate. (Set NT_GEN_LIB to a specific archive, or NT_ALLOW_STALE_LIB=1 to override — the latter is almost never right.)";

genLibStaleQ[lib_, incDir_] := Module[{hdrs, newest},
    hdrs = FileNames["*.hpp" | "*.h", incDir, Infinity];
    If[hdrs === {},
      Return[False]];
    newest = Max[AbsoluteTime /@ (FileDate[#, "Modification"]& /@ hdrs)];
    AbsoluteTime[FileDate[lib, "Modification"]] < newest];

resolveGenLib[incDir_] := Module[{env, base, cands, lib},
    env = Environment["NT_GEN_LIB"];
    lib =
      If[StringQ[env] && FileExistsQ[env],
        env,
        base = DirectoryName[incDir];(* <prefix> (installed) or <repo>/numtracer (in-tree): parent of include/ *)
        cands =
          {
            FileNameJoin[{base, "lib", "libNumTracer.a"}],
            (* installed layout *)
            FileNameJoin[{base, "build", "libNumTracer.a"}]
          };(* in-tree default build dir *)
        SelectFirst[cands, FileExistsQ, $Failed]];
    If[lib =!= $Failed && Environment["NT_ALLOW_STALE_LIB"] === $Failed && genLibStaleQ[lib, incDir],
      Module[{hdrs = FileNames["*.hpp" | "*.h", incDir, Infinity], newest},
        newest = First @ SortBy[hdrs, -AbsoluteTime[FileDate[#, "Modification"]]&];
        Message[MakeNTKernel::stalelib, lib, newest];
        Abort[]]];
    lib];

(* ---- semantic complexQ: probe whether the imaginary part actually vanishes -------------------
   The syntactic `complexQ = !FreeQ[Diagrams, Complex]` only sees that SOME diagram coefficient carries
   an `i` (projector i, non-abelian/quark colour f^abc T^bT^c = (iN/2)T^a). Those i's frequently pair up
   (i·i = -1) across diagrams so the assembled flow is exactly REAL — but the cancellation involves the
   concrete trace VALUES, and at the Mathematica stage the traces are opaque generated C++ symbols, so
   it can't be seen there. The just-generated C++ traces CAN see it: this compiles+runs a tiny probe that
   evaluates Im(integrand) over random frames with smooth real stub dressings/regulators. Returns True iff
   the imaginary part vanishes (|Im| <= tol·|Re|) over all sampled points — in which case the caller
   re-emits a real (double) kernel losslessly. Conservative: any failure (compile/run/no-points) -> False
   (keep the complex kernel). *)
(* "TraceArrayDecl": with CrossTraceCSE the integrand's trace tokens are `tarr[i]` reads, not
   `ns::tr_i(fenv)` calls, so the probe TU must declare and fill that array exactly as the kernel's
   coreBlock does — otherwise the probe fails to compile, the failure is (deliberately) swallowed as
   "keep the complex kernel", and the flow silently loses the lossless RePart double-kernel emission.
   Empty for the per-trace path. *)

Options[numericImagProbeRealQ] = {"NPoints" -> 4000, "Tol" -> 1.*^-9, "TraceArrayDecl" -> ""};

numericImagProbeRealQ[integrand_, args_, fillArgs_, angleDefs_, angleDecls_, nsHome_, headerFile_, drTable_ : <||>, opts : OptionsPattern[]] :=
  Module[{keepHeads, keepSyms, seedOf, argComb, stub, probeFull, probeProj, probeParams, probePre, fnFull, fnProj, drDecls, drFillArgs, randDecls, callArgs, src, cppFile, bin, rc, out, np, tol, tracesDir, cxx, parsed, distOf},
    np = OptionValue["NPoints"];
    tol = OptionValue["Tol"];
    tracesDir = DirectoryName[headerFile];
    cxx = resolveGenCxx[];
(* GENERAL stubbing: replace EVERY external real-valued atom — any dressing (any arity), any
   regulator/support function, any named constant — with `ntStub(seed_head, hash(args))`, an
   INDEPENDENTLY-SEEDED pseudo-random real per head. "External" = head is not a structural math
   operation, not a frame arg, not a known constant, and not a trace token (those are raw C++ strings,
   left untouched). Independence is essential: a single shared stub would make `dress1[x]-dress2[x]`
   collapse to 0 and mask a real surviving imaginary part. Distinct heads -> distinct seeds; distinct
   arguments -> distinct hashes -> distinct values, so any genuine imaginary part is exposed for
   ARBITRARY dressings. *)
    keepHeads = Alternatives[Plus, Times, Power, Rational, Sqrt, Sin, Cos, Tan, Cot, Sec, Csc, Exp, Log, Abs, Sign, ArcTan, ArcSin, ArcCos, Sinh, Cosh, Tanh, Max, Min, Floor, Ceiling, Mod, Complex, List, Global`ntStub];
    keepSyms = Join[args, First /@ angleDefs, {Pi, E, EulerGamma, Degree, GoldenRatio}];
    seedOf[h_] := N[Mod[Hash[SymbolName[h]], 100003] + 7];(* distinct seed per external head *)
    argComb[a_List] := Total[MapIndexed[#1 * N[GoldenRatio] ^ (#2[[1]] - 1)&, a]];(* 0 for a constant *)
    stub[e_] := Module[{x},
        x = e //. (h_Symbol)[a___] /; FreeQ[keepHeads, h] && !MemberQ[args, h] :> Global`ntStub[seedOf[h], argComb[{a}]];
        x /. (s_Symbol) /; !MemberQ[keepSyms, s] && FreeQ[keepHeads, s] && Context[s] =!= "System`" && s =!= Global`ntStub :> Global`ntStub[seedOf[s], 0.]
      ];
    probeFull = stub[integrand];
    probeProj = probeFull /. Complex[a_, b_] :> a;(* the candidate real projection *)
    probeParams = (<|"Name" -> SymbolName[#], "Type" -> "double", "Const" -> True, "Reference" -> True|>)& /@ args;
(* DRESSED kernels: the generated `fill()` takes one extra `double dr_<id>` per dressing atom (the
   kernel body computes each atom from regulators/interpolators and passes the VALUE). The probe has
   none of that runtime in scope, so it must compute each atom with the SAME pseudo-random stubbing
   used for the integrand and pass the values — otherwise the probe `.cpp` fails to compile (too few
   args to fill) and the kernel is conservatively kept COMPLEX, silently losing the lossless RePart
   double-kernel emission (the imaginary-half DCE that makes the distributed baseline fast). Atoms can
   reference the derived angles, so the angle decls precede them. *)
    drDecls =
      KeyValueMap[
        Function[{id, atom},
          "const double dr_" <> ToString[id] <> " = " <> cppFlat[stub[atom]] <> ";"],
        drTable];
    drFillArgs = ("dr_" <> ToString[#])& /@ Sort[Keys[drTable]];
    probePre =
      StringRiffle[
        Join[
          angleDecls,
          {"double fenv[(" <> nsHome <> "::nenv) > 0 ? (" <> nsHome <> "::nenv) : 1];"},
          drDecls,
          {nsHome <> "::fill(fenv, " <> StringRiffle[Join[SymbolName /@ fillArgs, drFillArgs], ", "] <> ");"},
          If[OptionValue["TraceArrayDecl"] === "",
            {},
            {OptionValue["TraceArrayDecl"]}]],
        "\n"];
    fnFull = FunKit`MakeCppFunction[probeFull, "Name" -> "probe_full", "Prefix" -> "static inline", "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> probeParams, "Body" -> probePre];
    fnProj = FunKit`MakeCppFunction[probeProj, "Name" -> "probe_proj", "Prefix" -> "static inline", "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> probeParams, "Body" -> probePre];
    distOf[a_] := Which[
        StringContainsQ[SymbolName[a], "cos"],
          "Uc",
        StringContainsQ[SymbolName[a], "phi"],
          "Uph",
        True,
          "U"];
    randDecls = StringRiffle[("double " <> SymbolName[#] <> " = " <> distOf[#] <> "(rng);")& /@ args, " "];
    callArgs = StringRiffle[SymbolName /@ args, ", "];
    src =
      StringJoin[
        "#define __host__\n#define __device__\n",
        "#include <complex>\n#include <cmath>\n#include <random>\n#include <cstdio>\n",
        "#include \"" <> FileNameTake[headerFile] <> "\"\n",
        "template<int N> static inline double powr(double x){ double r=1.0; int n=N<0?-N:N; for(int i=0;i<n;++i) r*=x; return N<0?1.0/r:r; }\n",
        "using std::pow; using std::sqrt; using std::sin; using std::cos; using std::tan; using std::exp; using std::log; using std::fma; using std::fabs;\n",
        "static inline std::complex<double> fma(const std::complex<double>&a,const std::complex<double>&b,const std::complex<double>&c){return a*b+c;}\n",
        "template<class T> using complex = std::complex<T>;\n",
(* independently-seeded pseudo-random real in [0.4,0.9): same (seed,arg) -> same value (a dressing is
   a function), distinct (seed,arg) -> independent value, so no two dressings or arguments collide. *)
        "static inline double ntStub(double seed, double x){ double h = std::sin(seed*0.1031 + x*0.3127 + 1.7)*43758.5453; return 0.4 + 0.5*(h - std::floor(h)); }\n",
        fnFull,
        "\n",
        fnProj,
        "\n",
        "int main(){ std::mt19937_64 rng(12345); std::uniform_real_distribution<double> U(0.25,3.0),Uc(-0.9,0.9),Uph(0.1,6.2);\n",
        "  double mim=0,mdiff=0,mre=0,mrim=0,mrdiff=0; long ok=0;\n",
        "  for(int n=0;n<" <> ToString[np] <> ";++n){ " <> randDecls <> "\n",
        "    std::complex<double> f = probe_full(" <> callArgs <> "); std::complex<double> pj = probe_proj(" <> callArgs <> ");\n",
        "    double im=std::imag(f), re=std::real(f), df=std::abs(f-pj);\n",
(* PER-POINT relative measures (mrim, mrdiff): a global max|Im|/max|Re| can let a localized
   imaginary part hide behind a large |Re| at some OTHER point (catastrophic cancellation). The
   +1 floor degrades gracefully to an absolute test when |Re|/|f| are small. The verdict keys on
   these; the absolute trio is kept only for the log. *)
        "    if(std::isfinite(im)&&std::isfinite(re)&&std::isfinite(df)){ mim=std::max(mim,std::fabs(im)); mdiff=std::max(mdiff,df); mre=std::max(mre,std::fabs(re));\n",
        "      mrim=std::max(mrim, std::fabs(im)/(std::fabs(re)+1.0)); mrdiff=std::max(mrdiff, df/(std::abs(f)+1.0)); ++ok; } }\n",
        "  std::printf(\"%.10e %.10e %.10e %.10e %.10e %ld\\n\", mim, mdiff, mre, mrim, mrdiff, ok); return 0; }\n"
      ];
    cppFile = FileNameJoin[{$TemporaryDirectory, "ntprobe_" <> StringReplace[nsHome, {":" -> "_"}] <> ".cpp"}];
    bin = StringReplace[cppFile, ".cpp" -> ""];
    ntExportCpp[cppFile, src];
    rc = Run[cxx <> " -std=c++20 -O1 -w -I '" <> tracesDir <> "' '" <> cppFile <> "' -o '" <> bin <> "' 2> '" <> bin <> ".cerr'"];
    If[rc =!= 0,
      ntLog["[probe] compile failed (rc=", rc, ") — keeping complex kernel (conservative)"];
      Return["Complex"]];
    Run["'" <> bin <> "' > '" <> bin <> ".out'"];
    out =
      If[FileExistsQ[bin <> ".out"],
        Import[bin <> ".out", "Text"],
        ""];
    parsed = Quiet @ Check[ToExpression[StringReplace[#, {"e+" -> "*^", "e-" -> "*^-", "e" -> "*^"}]]& /@ StringSplit[StringTrim[out]], $Failed];
    If[!MatchQ[parsed, {_?NumericQ, _?NumericQ, _?NumericQ, _?NumericQ, _?NumericQ, _?NumericQ}] || parsed[[6]] < 1,
      ntLog["[probe] run produced no usable points — keeping complex kernel"];
      Return["Complex"]];
    ntLog["[probe] over ", Round[parsed[[6]]], " pts:  max|Im|=", ScientificForm[parsed[[1]], 3], "  max|full-proj|=", ScientificForm[parsed[[2]], 3], "  max|Re|=", ScientificForm[parsed[[3]], 3], "  rel|Im|=", ScientificForm[parsed[[4]], 3], "  rel|full-proj|=", ScientificForm[parsed[[5]], 3]];
(* Three-way verdict from the C++ evaluation (which resolves every complex multiplication):
     "Complex"  Im survives                          -> genuinely complex, keep it.
     "Pure"     Im=0 AND Complex->Re projection exact -> drop imaginary coeffs (clean real arithmetic).
     "RePart"   Im=0 but projection differs           -> a trace is itself complex; the value is real
                                                         but only `.real()` of the full complex result
                                                         is correct, so wrap the return in a C++ real part. *)
    (* verdict on the PER-POINT relative measures (mrim, mrdiff): no global-scale inflation. *)
    Which[
      parsed[[4]] > tol,
        "Complex",
      parsed[[5]] <= tol,
        "Pure",
      True,
        "RePart"]];

(* ---- group-diagonal dressing fold: SUNPoly via the validated C++ engine ---------------------
   Each diag-dressed colour-net STRING (carrying sun<n>.diag{Fund,Adj}(...,{d0,…}) factors) is folded
   by sun_value_dressed in a tiny build-time program (the same emit/compile/run seam as the imaginary
   probe), returning per net a list of terms {coeffRe, coeffIm, {dr, ...}} (a flat list of dressing
   ids, repetition = power). Reuses the numeric engine verbatim, so the per-component colour weights
   are byte-identical to the typed-out SU(N) tables — no Mathematica reimplementation of the algebra. *)

diagColPolys[colnetStrs_, includeDir_] :=
  Module[{cxx = resolveGenCxx[], src, cppFile, bin, rc, out, lines, res = {}, cur = Null, num},
    num[s_] := ToExpression[StringReplace[s, {"e+" -> "*^", "e-" -> "*^-", "e" -> "*^"}]];
    src =
      StringJoin[
        "#include \"numtracer/network/sun_net.hpp\"\n#include <cstdio>\n#include <vector>\n",
        "using namespace numtracer; using namespace numtracer::network;\n",
        "int main(){\n",
        (* one SUNEnv per distinct rank in the diag colour nets (colFacG emits `sun<n>.diag…` factors). *)
        StringJoin["  SUNEnv sun" <> # <> "(" <> # <> ");\n"& /@ DeleteDuplicates @ Flatten @ StringCases[colnetStrs, "sun" ~~ r : DigitCharacter.. ~~ "." :> r]],
        "  std::vector<SUNNet> nets = {" <> StringRiffle[colnetStrs, ", "] <> "};\n",
        "  for(std::size_t n=0;n<nets.size();++n){\n",
        "    SUNPoly p = sun_value_dressed(nets[n]);\n",
        "    std::printf(\"NET %zu %zu\\n\", n, p.size());\n",
        "    for(const auto& t : p){\n",
        "      std::printf(\"T %.17g %.17g %zu\", t.coeff.re, t.coeff.im, t.dress.size());\n",
        "      for(int d : t.dress) std::printf(\" %d\", d);\n",
        "      std::printf(\"\\n\"); } }\n  return 0;\n}\n"];
    cppFile = FileNameJoin[{$TemporaryDirectory, "ntdiagpoly.cpp"}];
    bin = StringReplace[cppFile, ".cpp" -> ""];
    ntExportCpp[cppFile, src];
(* HEADER_ONLY: this one-off TU calls a SPLIT engine entry point (sun_value_dressed), whose body is
   `#if NUMTRACER_DEFINE_BODIES` — in a normal consumer TU only the declaration is visible and the
   definition is linked from libNumTracer.a. This helper links nothing, so without the define it
   fails at link time with an undefined reference (which aborted every diagonal-colour flow:
   gen_flavour_ingroup, gen_gluon_condensate). The TU is tiny, so inlining the bodies is free. *)
    rc = Run[cxx <> " -std=c++20 -O1 -w -DNUMTRACER_HEADER_ONLY=1 -I '" <> includeDir <> "' '" <> cppFile <> "' -o '" <> bin <> "' 2> '" <> bin <> ".cerr'"];
    If[rc =!= 0,
      Print["[diagpoly] compile failed (rc=", rc, "):\n", Quiet @ Check[Import[bin <> ".cerr", "Text"], ""]];
      Abort[]];
    Run["'" <> bin <> "' > '" <> bin <> ".out'"];
    out =
      If[FileExistsQ[bin <> ".out"],
        Import[bin <> ".out", "Text"],
        ""];
    lines = Select[StringSplit[StringTrim[out], "\n"], # =!= ""&];
    Do[
      Module[{tk = StringSplit[ln]},
        Which[
          tk[[1]] === "NET",
            If[cur =!= Null,
              AppendTo[res, cur]];
            cur = {},
          tk[[1]] === "T",
            Module[{re = num[tk[[2]]], im = num[tk[[3]]], m = ToExpression[tk[[4]]]},
              AppendTo[
                cur,
                {
                  re,
                  im,
                  If[m === 0,
                    {},
                    ToExpression /@ tk[[5 ;; 4 + m]]]}]]]],
      {ln, lines}];
    If[cur =!= Null,
      AppendTo[res, cur]];
    res];

(* The invariant-basis GENERATION path (the default numeric backend). Algorithm:
     1. reset the per-generation registries ($ctCache, resetDiagDr, resetDr) and unpack the NTKernel.
     2. per diagram, build the Lorentz/colour nets: splitColourGroupsInv groups branches by colour
        structure; a plain diagram lowers to a LorentzNet, a dressed-numerator one to a DPoly chain
        (the two paths the generator contracts differently). CSE-share repeated net/trace builders.
     3. emit a build-time C++ generator program (emitNumericGenerator) — a main TU plus per-net units —
        that runs reduce -> rebase -> lower at codegen time and PRINTS the straight-line traces header.
     4. compile + run that generator (resolveGenCxx), probe whether the result is real/complex to pick
        the kernel's number type, then emit and write (write-if-changed) the kernel header that fills
        the fundamental symbols and calls the generated trN(f). *)

mkGenerateKernel[NTKernel[k_], genFile_, kernelFile_, headerFile_, OptionsPattern[]] :=
  Module[{name, ns, dress, scalarParams, adParams, adNames, scalarTy, args, frame, env, nc, mask, ncomp, fillArgs, fillArgSig, invNets, invRest, g, colourNets, gcol, preamble, integrand, kernelParams, constParams, mkParam, kernelFn, constFn, classStr, header, hdrInc, incDir, genPre, genUnits, genDecl, genMain, declFile, unitFiles, genSrc, bin, run, hasFund, complexQ, colDecls, colToks, angleDefs, angleDecls, crossCSE, traceRef, nGrp, decor, tarrDecl, kns, sns, runInc, extraInc, interpTy, nsHome, gc, regTemplate, regAlias, symDefs = <||>, dmono = {}, atomStrs = {}, groupCombos = {}, groupContribs = {}, realOnlyG = {}, dressedIdx = {}, diagTokExpr = {}, factorNets = {}, lorFacOf = {}, pGroupOf = <||>, nAdd = 0, factorCompOf = <||>},
    Needs["FunKit`"];
(* A large flow assembles a kernel with one summand per diagram GROUP (ZA4: 1274). Several codegen
   steps (the integrand Sum, COEN's expression lowering) recurse ~linearly in that count, so the
   default $RecursionLimit of 1024 is exceeded — and $RecursionLimit::reclim does NOT abort, it
   returns a held expression and the script continues to "DONE" having SILENTLY skipped the kernel
   (the ZA4 silent-skip). Raise the limit generously for the emission; a real runaway would still
   hit the (much higher) ceiling and surface. *)
    $RecursionLimit = Max[$RecursionLimit, 1048576];
    name = OptionValue["Name"];
    dress = OptionValue["Dressings"];
    scalarParams = OptionValue["ScalarParams"];(* loop-independent scalar doubles threaded into the signature *)
(* AD-flagged scalars (d1V, d2V for FE-potential flows) must be `const auto&` so the kernel also
   accepts autodiff::real from the integrator_AD twin; everything else stays `const double&`. *)
    adParams = OptionValue["ADParams"];
    adNames =
      If[StringQ[#],
          #,
          SymbolName[#]
        ]& /@ adParams;
    scalarTy =
      Function[nm,
        If[MemberQ[
          adNames,
          If[StringQ[nm],
            nm,
            SymbolName[nm]]],
          "auto",
          "double"]];
    angleDefs = OptionValue["AngleDefs"];
    crossCSE = OptionValue["CrossTraceCSE"];
    gc = OptionValue["GlobalCollect"];                                                              (* Route-B: fold the WHOLE integrand (all diagrams, dressings as
inert symbols) into ONE collected polynomial -> one kernel, like FORM. crossCSE is subsumed by it. *)
    decor = OptionValue["Decorator"];
    kns = OptionValue["KernelNamespace"];
    sns = OptionValue["SupportNamespace"];
    runInc = OptionValue["RuntimeInclude"];
    extraInc = OptionValue["ExtraIncludes"];
(* the fRG/DiFfRG kernel shape (template<typename REG> + the private REG:: wrappers) is opt-in; the
   Regulator alias is meaningless without it, so it implies it. *)
    regAlias = TrueQ[OptionValue["RegulatorAlias"]];
    regTemplate = TrueQ[OptionValue["RegulatorTemplate"]] || regAlias;
    interpTy = ntDressType[OptionValue["DressingType"]];
    ns = OptionValue["Namespace"] /. Automatic -> ToLowerCase[name];
    nsHome = kns <> "::" <> ns;(* where the generated trace fns / nenv / fill live *)
    args = k["Args"];
    frame = k["Frame"];
    env = k["Env"];
(* Vestigial compile-chain salt: every SU(N) head carries its own rank N, so builders read it from
   the head and never consult this positional arg (kept only to avoid re-threading the recursion). *)
    nc = 0;
    mask = Association @ KeyValueMap[#1 -> frameMask[resolveComponents[#1, frame]]&, env];
    fillArgs = Select[args, # =!= Global`k&];(* scalars the fill needs *)
(* NUMERIC backend: build the component table over a compact parametrisation. Automatic uses the
   unit-loop spec (loop = magnitude × unit-direction symbols + unit constraint) so the contraction is
   as compact as the sp basis for BOTH Lorentz and Dirac nets. User "Components" are taken verbatim
   (polynomialised for any Sqrt/trig). *)
    ncomp =
      Module[{uc = OptionValue["Components"], ud = OptionValue["SymbolDefs"], pf, ad, ug},
        {pf, ad, ug} =
          Which[
            uc =!= Automatic,
              Append[polyFrameSpec[uc], {}],
            (* explicit user Components *)unitLoopOkQ[frame, Global`p, Global`l1],
              unitLoopFrameSpec[frame, Global`p, Global`l1],
            (* compact vacuum SP frame *)True,
              Append[polyFrameSpec[frame], {}]
          ];(* finite-T / general frame *)
        symDefs = Join[ad, ud];
        numericComponents[env, pf, symDefs, ug]];
(* [[maybe_unused]]: a frame may not reference every fill() argument (e.g. an angle or dressing atom
   that only some diagrams use), so mark each parameter to keep the emitted kernel -Wunused-clean. *)
    fillArgSig = StringRiffle[("[[maybe_unused]] double " <> SymbolName[#])& /@ fillArgs, ", "];
(* walk diagrams: each one is a Lorentz trace x a colour factor x a dressing/kinematic coeff. The
   bridge distribution split each 4-gluon vertex into colour channels, so MANY diagrams share the
   SAME coeff and differ only in (colour x Lorentz). Collect per diagram {coeff, lorNet, colNet},
   then GROUP by coeff: per group the generator folds the (constant) colour into the Lorentz poly
   and combines them — combinedTr_g = Σ_d colN_d · trN_d — so the kernel evaluates ~one polynomial
   per Feynman graph (≈5), not one per channel (51). Vanishing-colour diagrams drop out for free. *)
(* Fundamental colour (a quark loop's T^a) can't go through the adjoint SUNNet, but it IS a
   small et-instantiable contraction (unlike the four-gluon type), so we evaluate it with the et
   engine as a kernel `constexpr double` — exactly the direct path's constant-colour handling —
   and multiply it into the integrand. The generator then carries an identity colour (SUNNet{}).
   With fundamental colour present we do NOT group diagrams (each keeps its own colour constant). *)
    hasFund = !FreeQ[k["Diagrams"], _ntSUNT | _ntSUNDeltaFund];
(* Dirac VERTEX: projector i + imaginary non-abelian colour (f^abc T^b T^c = (iN/2)T^a). Per
   diagram coeff*colour is real (the i's combine), so we keep the colour constant COMPLEX and
   take Re of the assembled integrand. *)
    complexQ = !FreeQ[k["Diagrams"], Complex];
(* CrossTraceCSE (IMPLEMENTED 2026-07-19; it was a documented stub before — `trace_all` existed only
   as the call site emitted below and in comments, so turning the option on produced a kernel that
   failed to compile with "'trace_all' is not a member of ...", and the complexQ guard that used to
   sit here was masking that).

   What it targets: cross-trace sharing is the one mechanism that attacks the measured redundancy in
   these flows — ZAqbq1_147 Mq-in carries 28,856 monomials across its traces of which only 8,144 are
   DISTINCT (3.54x; some recur in 24 different traces), and FormTracer's advantage here is exactly
   that it sums all diagrams into ONE polynomial before expanding, so those duplicates collect
   (measured Mq-out: FORM 5,060 monomials vs NumTracer 8,482, BOTH lowering at an identical 2.14
   ops/monomial). The compiler cannot do it for us: collecting one monomial out of N traces is a
   floating-point reassociation across function boundaries (force-inlining all 108 traces recovered
   only -3.5% instructions). See example/ZAqbq147-MqBench/FINDINGS.md.

   The old complex restriction is GONE. `tarr` is emitted as std::complex<double> whenever any trace
   is complex — the same type trN(fenv) returned — so the kernel's ntRe/ntIm reads are unchanged and
   nothing is truncated. Do NOT reintroduce a real tarr with per-trace phase tracking. *)
    invNets = {};
    invRest = {};
    colourNets = {};
    colDecls = {};
    colToks = {};
    preamble = {};
    factorNets = {};
    lorFacOf = {};
    factorCompOf = <||>;(* factor net indices, per-net factor-id list, net->factor-id *)
    $ctCache = <||>;(* clear the compileTInv memo cache for this generation *)
(* dressedSlotStr / orderDiracLoops memos: both depend on generation-fixed state ($ntDressResolve,
   env), so they MUST be cleared here alongside $ctCache. *)
    $dsCache = <||>;
    $odCache = <||>;
    resetDiagDr[];(* clear the per-component diagonal-dressing registry for this generation *)
    resetDr[];(* clear the scalar-dressing (ntDressedNum) registry for this generation *)
(* frame resolver for dressed-numerator option coefficients (compileDirac → dressedSlotStr): the same
   ntSP/ntSPS/ntVec[q,i] → component substitution used for diag["Coeff"] below. *)
    $ntDressResolve =
      Function[s,
        s /. {ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame], ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]], ntVec[q_, ii_Integer] :> resolveComponents[q, frame][[ii + 1]]}
      ];
    Module[{diagData = {}},
      ntLog[
        "[prof] per-diagram net-build (",
        Length[k["Diagrams"]],
        " diagrams): ",
        First @
          AbsoluteTiming[
            MapIndexed[
              Function[{diag, di},
                Module[{coeff, colBr, constAcc = {}, colTok = "", entries = {}, d = di[[1]] - 1, pureLorAcc = {}, nNonConst = 0, nNCDirCol = 0, diracComps = {}},
                  coeff = diag["Coeff"] /. {ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame], ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]], ntVec[q_, i_Integer] :> resolveComponents[q, frame][[i + 1]]};
                  MapIndexed[
                    Function[{comp, ci},
                      Module[{},
                        If[comp["Constant"],
(* Constant SU(N) component — colour and/or flavour, fundamental and/or adjoint: every group
   head carries its own rank N, so all fold through the SINGLE numeric SUNNet path
   (sun_value_cx, contracting each rank separately and multiplying). The fold is COMPLEX
   (-> Cx): when an imaginary non-abelian-vertex colour (f^abc T^b T^c = (iN/2) T^a) is folded
   in, the diagram's trace lowers to a complex `tr_i`; the kernel multiplies by the complex
   dressing coefficient and the consumer takes Re, so the imaginary part survives. ACCUMULATE
   components — a Yukawa loop carries SEVERAL constant components (a colour trace AND a flavour
   trace); mergeColNet folds them together so none is dropped (a bare `col = str` would keep
   only the last, e.g. the δ^ii = Nf factor, giving a ~50% wrong trace).
   ACCUMULATE THE FACTORS and compile once after the loop (see colBr below) rather than
   compiling per component and splicing the strings with mergeColNet: a constant component
   may be a PLUS (the four-quark Fierz flavour structure δδ - 4·T·T is one), which has no
   single-net representation, and one Expand of the whole constant product does the
   cross-product across several summed components in one step. Compiling the product is
   equivalent to the old per-component splice — mergeColNet only concatenates factor
   lists — so single-branch flows regenerate byte-identically. *)
                          constAcc = Join[constAcc, comp["Factors"]],
(* Non-constant component. The DISCONNECTED components of ONE diagram MULTIPLY (their scalar
   trace values `Times @@ toks`) — they are NOT separate summed diagrams. Collect them so
   the post-loop assembly can form the
   product (see there). Route by structure:
     - any colour (entangled in a Plus, or a top-level T^a × …) or a gamma chain: collect the
       component's splitColourGroupsInv entries (dirac_value net per colour branch / chunked
       Lorentz net) as ONE Dirac/colour component in diracComps;
     - pure Lorentz (no colour, no gamma): accumulate factors — ALL pure-Lorentz components
       fold into ONE product net (disjoint ids make the C++ contract_factors multiply them). *)
                          (
                            nNonConst++;
                            If[colourEntangledQ[comp["Factors"]] || !FreeQ[comp["Factors"], _ntGamma | _ntGamma5 | _ntDeltaDirac | _ntDressedNum | _ntDiracSlot],
                              (
                                nNCDirCol++;
                                AppendTo[diracComps, splitColourGroupsInv[comp["Factors"], diag["Ids"], env, mask, nc]]
                              ),
                              pureLorAcc = Join[pureLorAcc, comp["Factors"]]])]]],
                    diag["Components"]];
(* The diagram's CONSTANT colour/flavour part, as a list of {netString, scalar} branches — one
   branch unless a constant component was a sum. Colour folds to a scalar (sun_value_cx -> Cx)
   and the generator already sums colour by emitting several nets into one group, so a summed
   colour component costs one extra net record per branch and needs no C++ support. *)
                  colBr =
                    If[constAcc === {},
                      {{"SUNNet{}", 1}},
                      compileColGSum[Times @@ constAcc, diag["Ids"]]];
(* ---- assemble the diagram's nets from its non-constant components -----------------------------
   A diagram with K disconnected non-constant components is a PRODUCT of K independent closed
   scalars (each a Dirac/colour trace or a pure-Lorentz scalar): `coeff * Times @@ toks`.
   ALL pure-Lorentz components fold into ONE product factor; each
   Dirac/colour component is its own factor. Two regimes:
     - <= 1 non-constant factor: the EXISTING additive path. The single Dirac component's entries
       (colour folded, GlobalCollect-fusible) OR the single combined pure-Lorentz product net is
       appended with diagData = coeff*scal. Flows without a disconnected diagram regenerate
       BYTE-IDENTICAL (any K>=2 diagram previously aborted, so none is committed).
     - >= 2 factors: FACTORED product. One Dirac component is the additive BASE (traceRef[gi],
       diagData = coeff*scal); every other component is emitted as its own fused trace GROUP
       (its per-entry scalar folded into the net, its colour folded by the group sum) and tagged
       with a factor id. The base carries the list of factor ids; the assembly multiplies in
       Π traceRef[factor groups] — P computed ONCE per component, no trace-polynomial blow-up.
   scalarleakCheck: each entry's restNet scalars (e[[4]] {restStr,scal}) become the generator's
   `dsc[]` numeric constants — a symbolic Lorentz tensor the net builder failed to fold would be
   CForm'd into undeclared C++. Catch it loudly, with the offender. *)
                  Module[{nDir = Length[diracComps], hasLor = pureLorAcc =!= {}, factorComps, factorIds = {}, scalarleakCheck, appendRec},
                    scalarleakCheck =
                      Function[es,
                        Do[
                          Module[{badS = FirstCase[ee[[4]], {_, s_} /; !NumericQ[s] :> s, Missing[]]},
                            If[!MissingQ[badS],
                              Message[mkGenerateKernel::scalarleak, d, badS];
                              Abort[]]],
                          {ee, es}]];
                    (* rec = {colNet, cores, dData, restList, lorFac(None|{ids..}), factorId(None|id)}; returns net idx *)
                    appendRec =
                      Function[rec,
                        Module[{ni = Length[invNets]},
                          invNets = Append[invNets, rec[[2]]];
                          invRest = Append[invRest, rec[[4]]];
                          colourNets = Append[colourNets, rec[[1]]];
                          colToks = Append[colToks, colTok];
                          diagData = Append[diagData, rec[[3]]];
                          lorFacOf = Append[lorFacOf, rec[[5]]];
                          If[rec[[6]] =!= None,
                            factorNets = Append[factorNets, ni];
                            factorCompOf[ni] = rec[[6]]];
                          ni]];
                    If[nDir + Boole[hasLor] <= 1,
                      (* ---- single non-constant factor (or none): EXISTING additive path, byte-identical ---- *)
                      Module[{
                        baseEntries =
                          Which[
                            nDir == 1,
                              diracComps[[1]],
                            hasLor,
                              ({"SUNNet{}", {#[[1]]}, #[[2]], {{"", 1}}}& /@ chunkLorInv[Times @@ pureLorAcc, diag["Ids"], env, mask, nc]),
                            True,
                              {}]},
                        scalarleakCheck[baseEntries];
                        Do[
                          appendRec[{mergeColNet[cb[[1]], e[[1]]], e[[2]], coeff cb[[2]] e[[3]], e[[4]], None, None}],
                          {cb, colBr},
                          {
                            e,
                            If[baseEntries === {},
                              {{"SUNNet{}", {"konst(1.0)"}, 1, {{"", 1}}}},
                              baseEntries]}]],
(* ---- >= 2 factors: factored product of disconnected components ----
   EVERY non-constant component becomes its OWN fused trace group: its colour-branch entries
   are summed WITHIN the group (GlobalCollect-style, colour folded by the group sum, the entry
   scalar folded into the net) so the group trace IS the component scalar. The diagram then
   contributes ONE additive ANCHOR term = coeff * colv(col) * Π(component group traces).
   Crucially this does NOT force the components' many entries into singletons (which would
   defeat colour-channel fusion and explode the trace count for a four-quark-dressed loop) —
   each component is exactly ONE trace (one scalar per component). *)
                      (
                        factorComps =
                          If[hasLor,
                            Append[diracComps, {{"SUNNet{}", {#[[1]]}, #[[2]], {{"", 1}}}}&[compileTInv[Times @@ pureLorAcc, diag["Ids"], env, mask, nc]]],
                            diracComps];
                        Do[
                          Module[{compEntries = factorComps[[ci]], fid = Length[invNets]},
                            scalarleakCheck[compEntries];
                            AppendTo[factorIds, fid];
                            Do[appendRec[{e[[1]], e[[2]], 1, ({#[[1]], #[[2]] e[[3]]}&) /@ e[[4]], None, fid}], {e, compEntries}]
                          ],
                          {ci, 1, Length[factorComps]}];
(* the anchor: a trivial unit net carrying the diagram coeff, the constant colour branch
   (folded once, via the group sum colv(net)*1), and the list of component factor ids.
   ONE anchor per colour branch — each is an independent additive term
   coeff·scal_b·colv(net_b)·Π(factor traces), and they all reference the SAME factorIds,
   so the expensive component traces are computed once and shared. *)
                        Do[appendRec[{cb[[1]], {"konst(1.0)"}, coeff cb[[2]], {{"", 1}}, factorIds, None}], {cb, colBr}]
                      )]]]],
              k["Diagrams"]]],
        " s"];
(* memo sizes: the net-build's cost is dominated by the DISTINCT Dirac/Lorentz structures, not the
   call count (a dense flow calls these tens of thousands of times for a few hundred distinct
   results). If a flow ever shows these growing with the call count, a memo has stopped hitting. *)
      ntLog["[prof]   memo sizes: compileTInv ", Length[$ctCache], " | orderDiracLoops ", Length[$odCache], " | dressedSlotStr ", Length[$dsCache]];
(* ---- per-component diagonal dressings (ntSUNDiag{Fund,Adj}) ----------------------------------
   A diagram whose colour net carries a diag factor folds (via the validated C++ engine,
   sun_value_dressed, run through the build-time seam) to a SUNPoly Σ_t coeff_t Π D^{dr}, where
   each dr names a distinct SCALAR runtime dressing (the surviving/named components; dropped ones
   vanish from the sum). Since the dressings are runtime, the colour can't be a compile-time
   constant baked into the trace: instead we (a) replace the diagram's colour net by the IDENTITY
   (so the generator folds colv=1 and the trace stays colour-free), and (b) build a runtime token
   `Σ_t coeff_t Π name(scale)` — ordinary scalar-dressing tokens — multiplied into the integrand.
   The Dirac/Lorentz trace is still computed ONCE, not a diagram per component. *)
      diagTokExpr = Table[1, {Length[colourNets]}];
      dressedIdx = Select[Range[Length[colourNets]], StringContainsQ[colourNets[[#]], ".diag"]&];
      If[dressedIdx =!= {},
        Module[{polys, resolveScale, incDir = OptionValue["IncludeDir"] /. Automatic :> resolveIncludeDir[]},
          resolveScale[s_] := s /. {ntSP[x_, y_] :> resolveComponents[x, frame] . resolveComponents[y, frame], ntSPS[x_, y_] :> Rest[resolveComponents[x, frame]] . Rest[resolveComponents[y, frame]], ntVec[q_, ii_Integer] :> resolveComponents[q, frame][[ii + 1]]};
          polys = diagColPolys[colourNets[[dressedIdx]], incDir];
          MapThread[
            Function[{d, p},
              diagTokExpr[[d]] =
                Total[
                  Function[term,
                      (term[[1]] + I term[[2]]) *
                        (
                          Times @@
                            (
                              Function[dr,
                                  $diagDrTable[dr]["Name"][resolveScale[$diagDrTable[dr]["Scale"]]]
                                ] /@ term[[3]]))
                    ] /@ p];
              colourNets[[d]] = "SUNNet{}"],
            {dressedIdx, polys}]];
        ntLog["[prof] diagonal-dressed diagrams: ", Length[dressedIdx], " (per-component colour-sum folded via sun_value_dressed seam)"]
      ];
(* trace reference: with CrossTraceCSE the kernel fills a `tarr[]` once via trace_all() and reads
   tarr[i]; otherwise it calls the independent tr_i(fenv). *)
      traceRef =
        If[crossCSE,
          "tarr[" <> ToString[#] <> "]",
          nsHome <> "::tr" <> ToString[#] <> "(fenv)"
        ]&;
(* Group diagrams (0-based) into traces by colour kind. ADJOINT-colour diagrams (colTok=="" —
   their real colour was folded numerically into the generator polynomial) FUSE by identical
   dressing coeff, so the kernel evaluates ~one polynomial per Feynman graph; cross-diagram like
   terms then collect in the invariant basis (many channels collapse to a few traces this way).
   FUNDAMENTAL-colour diagrams (a complex generator T^a, kept in the kernel as _colN so the trace
   polynomial stays REAL) stay singleton, each times its own _colN. A flow that mixes both colour
   kinds thus fuses the adjoint part and keeps the fundamental part per-diagram. *)
(* diag-dressed diagrams (diagTokExpr =!= 1) carry a per-diagram RUNTIME colour-sum token, so —
   like the fundamental-colTok diagrams — they stay singletons (they cannot fuse by dressing
   coefficient alone, the colour token differs). lorFac diagrams (lorFacOf =!= None — a Dirac trace
   times a disconnected pure-Lorentz scalar) likewise stay singletons: each carries a per-diagram
   multiplicative trace, so it must not fuse with another diagram's entries.
   The FACTOR nets (P, indices in factorNets) are EXCLUDED from the additive groups and appended as
   their own singleton trace groups at the tail of g — generated as traces but referenced only
   multiplicatively via lorFac, never summed into the integrand. nAdd marks the additive/factor
   boundary; pGroupOf maps a P net (0-based) to its (0-based) trace-group ordinal. *)
      Module[{adj, fund, additivePos, factorPos = (# + 1)& /@ factorNets, gAdd, gFactor},
        additivePos = Complement[Range[Length[diagData]], factorPos];
        adj = Select[additivePos, colToks[[#]] === "" && diagTokExpr[[#]] === 1 && lorFacOf[[#]] === None&];
        fund = Select[additivePos, colToks[[#]] =!= "" || diagTokExpr[[#]] =!= 1 || lorFacOf[[#]] =!= None&];
        gAdd = Join[(# - 1)& /@ GatherBy[adj, diagData[[#]]&], List /@ (fund - 1)];
(* each disconnected factor COMPONENT (possibly several colour-branch nets) fuses into ONE trace
   group, so its group trace = the component scalar (colour folded by the group sum). *)
        gFactor = GatherBy[factorNets, factorCompOf[#]&];(* factorNets are 0-based net indices *)
        g = Join[gAdd, gFactor];
        nAdd = Length[gAdd];
        pGroupOf = Association[MapIndexed[(factorCompOf[#1[[1]]] -> (nAdd + #2[[1]] - 1))&, gFactor]]];
(* GlobalCollect folds ALL colour numerically (colToks all empty above), so this groups EVERY
   diagram by its dressing coefficient — the quark-loop colour channels that the legacy path kept
   as fundamental singletons now MERGE into their Feynman graph, collapsing the trace count toward
   FORM's handful. The dressing coefficient stays FACTORED in `diagData` (COEN CSEs it, like FORM's
   _repl), so each group is one collected kinematic trace × its dressing — not a flat polynomial. *)
(* Sum the ADDITIVE groups only (1..nAdd); factor groups (nAdd+1..) are referenced multiplicatively
   via lorFac. A factored (disconnected) diagram multiplies in its other components' scalars
   Π traceRef[factor groups] (each computed ONCE, as a separate trace): the diagram's
   `coeff * Times @@ component-scalars`. *)
      integrand =
        Sum[
          With[{rep = g[[gi, 1]]},
            diagData[[rep + 1]] *
              If[colToks[[rep + 1]] === "",
                1,
                colToks[[rep + 1]]
              ] * diagTokExpr[[rep + 1]] *
              If[lorFacOf[[rep + 1]] === None,
                1,
                Times @@ (traceRef[pGroupOf[#]]& /@ lorFacOf[[rep + 1]])
              ] * traceRef[gi - 1]],
          {gi, nAdd}]];
(* OPTIMIZATION (opt-in, "PruneRealTraces"): a group whose dressing coefficient is REAL has only
   Re(trace) consumed (the consumer takes Re of the whole kernel; a real coeff cannot move
   Im(trace) into the real part). Flag it so the generator emits a `double` trace and never
   computes the dead imaginary half. Empty list (default) => generator emits all-complex traces,
   which the probe needs to verify cancellation — see the "PruneRealTraces" note. *)
    realOnlyG =
      If[TrueQ[OptionValue["PruneRealTraces"]],
        (FreeQ[diagData[[#[[1]] + 1]], Complex])& /@ g,
        {}];
    nGrp = Length[g];
(* The tarr declaration+fill, used by BOTH the kernel's coreBlock and the RealProbe TU (which
   evaluates the same integrand, so it needs the same tokens in scope). `trace_all_t` is emitted by
   emit_cpp_fused from the ACTUAL lowered roots (complex iff some trace is complex), so the array
   type can never disagree with what trace_all stores — and ntIm(double)=0.0 is then correct rather
   than lossy, because the type is double only when every trace really is real. *)
    tarrDecl = nsHome <> "::trace_all_t tarr[" <> ToString[nGrp] <> "]; " <> nsHome <> "::trace_all(fenv, tarr);";
(* [B2 localize] surface the post-net-build shape so a silent no-output (e.g. empty nets / empty
   grouping) is visible rather than appearing as a clean DONE. *)
    ntLog["[prof] post-net-build: nets=", Length[invNets], " groups(nGrp)=", nGrp, " complexQ=", complexQ];
    If[Length[invNets] === 0 || nGrp === 0,
      Message[mkGenerateKernel::emptynets, name, Length[invNets], nGrp];
      Abort[]];
    (* kinematic angle defs (kept symbolic in the dressing): emit once as named temporaries. *)
    angleDecls = ("const double " <> SymbolName[First[#]] <> " = " <> cppFlat[Last[#]] <> ";")& /@ angleDefs;
(* NB: deliberately NO `using std::complex;` — unqualified complex<double> resolves to the
   support namespace's `complex` alias. That indirection lets a device/CUDA support header
   substitute a device-safe complex (std::complex arithmetic lowers to gcc _Complex builtins
   that nvcc silently miscompiles to 0 in device code), without changing the emitted kernel. *)
(* the fenv setup block: declare fenv, (dressed only) compute each dressing atom into dr_<id>, fill,
   and (CrossTraceCSE) precompute the traces. *)
    With[{hasDr = !FreeQ[invNets, _ntDressedCore]},
      Module[{coreBlock},
        coreBlock =
          {
            "double fenv[(" <> nsHome <> "::nenv) > 0 ? (" <> nsHome <> "::nenv) : 1];",
            Sequence @@
              If[hasDr,
                KeyValueMap[
                  Function[{id, atom},
                    "const double dr_" <> ToString[id] <> " = " <> cppFlat[atom] <> ";"],
                  $drTable],
                {}],
            With[{
              fillCallArgs =
                If[hasDr,
                  Join[SymbolName /@ fillArgs, ("dr_" <> ToString[#])& /@ Sort[Keys[$drTable]]],
                  SymbolName /@ fillArgs]},
              nsHome <> "::fill(fenv, " <> StringRiffle[fillCallArgs, ", "] <> ");"],
            If[crossCSE,
              tarrDecl,
              Nothing]};
(* DRESSED: the dr_<id> dressing expressions can reference the derived kinematic angles, so the
   angle (and colour) decls must precede the fenv block. NON-dressed: keep the original order
   (fenv before angle/colour decls) so those kernels regenerate byte-identical. *)
        preamble =
          StringRiffle[
            If[hasDr,
              Join[ntSupportUsings[sns], angleDecls, colDecls, coreBlock, preamble],
              Join[ntSupportUsings[sns], coreBlock, angleDecls, colDecls, preamble]],
            "\n"]]];
    mkParam[nm_, ty_] := <|
        "Name" ->
          If[StringQ[nm],
            nm,
            SymbolName[nm]],
        "Type" -> ty,
        "Const" -> True,
        "Reference" -> True
      |>;
(* every dressing — including the named per-component diagonal dressings (ntSUNDiag{Fund,Adj}) —
   is an ordinary scalar interpolator kernel parameter. *)
    With[{
      dressTy =
        Function[nm,
          interpTy]},
      kernelParams = Join[mkParam[#, "double"]& /@ args, mkParam[#, scalarTy[#]]& /@ scalarParams, mkParam[#, dressTy[#]]& /@ dress];
      constParams = Join[mkParam[#, "double"]& /@ Select[args, # === Global`p || # === Global`k&], mkParam[#, scalarTy[#]]& /@ scalarParams, mkParam[#, dressTy[#]]& /@ dress];
(* dressed kernels: fill() takes one `double dr_<id>` per dressing atom — the kernel body computes
   the atom's value (regulators / interpolators in scope there) and passes it. Matches fm.dress. *)
      If[!FreeQ[invNets, _ntDressedCore],
        fillArgSig = fillArgSig <> StringJoin[(", [[maybe_unused]] double dr_" <> ToString[#])& /@ Sort[Keys[$drTable]]]
      ]];
(* LOUD GUARD: the integrand must be numeric-valued before it is lowered to C++. A DEGENERATE input
   — most often a basis whose Gram is singular at the chosen kinematics, so its inverse metric (and
   hence every dual projector) carries 0/0 — leaves Indeterminate / ComplexInfinity / DirectedInfinity
   in the coefficients. FunKit's lowering happily prints those as bare identifiers, emitting e.g.
   `return Indeterminate;` into the kernel header. That is a Mathematica symbol in generated C++:
   here it fails to compile only by luck (no such identifier), and a differently-named leak could
   compile into a silently wrong kernel.
   Observed with the FULL AqbqDirect basis (12 structures) at the symmetric point: Det[g] = 0 there,
   because the 12 structures are linearly DEPENDENT at that kinematic configuration — which is why
   the flows project with a restricted sub-basis. Refuse rather than emit. *)
    With[{bad = Cases[integrand, Indeterminate | _DirectedInfinity | ComplexInfinity, {0, Infinity}]},
      If[bad =!= {},
        Print["[NumTracer] ERROR: the integrand is not numeric — it contains ", Length[bad], " ", "Indeterminate/Infinity value(s), which would be emitted as bare Mathematica symbols in ", "the generated C++ (e.g. `return Indeterminate;`).\n", "  This almost always means a SINGULAR Gram at the chosen kinematics: the basis's ", "structures are linearly dependent there, so the inverse metric — and every dual ", "projector built from it — carries 0/0. Check Det[TBGetMetric[basis]] under the frame's ", "kinematics (e.g. the symmetric point), and project with a restricted sub-basis whose ", "Gram is non-degenerate. Offending value(s): ", Short[DeleteDuplicates[bad], 4]];
        Abort[]]];
(* the integrand -> C++ lowering (FunKit). Timed separately: it is the one heavy stage between the
   net-build and the generator emit, so without this the [prof] trail has a blind spot. *)
    ntLog[
      "[prof] FunKit kernel/class/header lowering: ",
      First @
        AbsoluteTiming[
          kernelFn = FunKit`MakeCppFunction[integrand, "Name" -> "kernel", "Prefix" -> decor, "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> kernelParams, "Body" -> preamble];
          constFn = ntConstFn[OptionValue["Constant"], decor, constParams, sns];
          classStr = ntKernelClass[name, {kernelFn, constFn}, decor, regTemplate, regAlias, {}];
          hdrInc = FileNameTake[headerFile];
          header =
            FunKit`MakeCppHeader[
(* the numeric kernel is flat straight-line arithmetic: the generated trace functions (hdrInc) plus
   the support runtime; no tensor-engine headers. *)"Includes" -> Join[extraInc, ntRuntimeIncludes[runInc], {"numtracer/sun/sun_data.hpp", hdrInc}], "Body" -> ntWrapBody[kns, classStr, name]
            ];],
      " s"];
    (* emit the generator source (the numeric matrix-product backend is the single generation path). *)
    ntLog["[prof] emitNumericGenerator: ", First @ AbsoluteTiming[{genPre, genUnits, genDecl, genMain} = emitNumericGenerator[invNets, invRest, colourNets, g, ncomp, ns, fillArgSig, kns, complexQ, realOnlyG, crossCSE];], " s"];
(* Split generator: a main TU + N net-builder unit TUs + a decl header (all in the tests/gen/ dir), so the
   net-builder codegen compiles in parallel (see emitNumericGenerator). The main `#include`s the decl. *)
    declFile = StringReplace[genFile, ".cpp" -> "_nets.hh"];
    unitFiles = Table[StringReplace[genFile, ".cpp" -> "_u" <> ToString[u - 1] <> ".cpp"], {u, 1, Length[genUnits]}];
    ntLog[
      "[prof] write generator files (",
      Length[unitFiles] + 2,
      " files, ",
      Round[(Total[StringLength /@ genUnits] + StringLength[genDecl] + StringLength[genMain]) / 1000000.],
      " MB): ",
      First @
        AbsoluteTiming[
          ntExportCpp[declFile, genDecl];
(* each unit #includes the shared decl header so its net builders can call the cross-unit CSE
   accessors (lc<k>()/dc<k>()) and sibling net builders, parsed once per TU. *)
          Module[{uInc = "#include \"" <> FileNameTake[declFile] <> "\"\n"},
            Do[ntExportCpp[unitFiles[[u]], uInc <> genUnits[[u]]], {u, 1, Length[genUnits]}]];
          genSrc = genPre <> "\n#include \"" <> FileNameTake[declFile] <> "\"\n\n" <> genMain;
          ntExportCpp[genFile, genSrc];],
      " s"];
    Print["wrote generator: ", genFile, " (+ ", Length[genUnits], " net units + decl header)"];
(* run the generator at codegen time -> the committed straight-line kernel header. The binary's
   stdout is redirected straight to the header FILE via the shell (Run), not captured in memory
   by RunProcess — the A4 kernel header is ~40k lines, and in-memory capture is fragile at that
   size. *)
    If[OptionValue["RunGenerator"],
      incDir = OptionValue["IncludeDir"] /. Automatic :> resolveIncludeDir[];
      bin = FileNameJoin[{$TemporaryDirectory, "gen_" <> ns}];
(* Time COMPILE and RUN separately — both count toward generation time, but the levers differ
   (compile: TU size / templates; run: reduce+rebase). Reported so neither is hidden. Compile the
   main TU (-O2) and the net-builder units (-O0) CONCURRENTLY (`&` + `wait`), then link — the unit
   codegen runs across cores. A failed unit compile leaves its .o missing -> the link rc is nonzero,
   so the existing rc check catches it. *)
      Module[
        {tcc, cc, mainObj, unitObjs, pcmd, lcmd, clog = bin <> "_compile.log", cxx = resolveGenCxx[], mainOpt, libPath = resolveGenLib[incDir], useLib, hoDef, libArg},
(* Default: link the prebuilt libNumTracer.a (engine bodies compiled once). If it is not found,
   fall back to a self-contained header-only compile so generation still works (older, slower
   path — every engine body re-instantiated in the main TU). *)
        useLib = StringQ[libPath] && FileExistsQ[libPath];
        hoDef =
          If[useLib,
            " ",
            " -DNUMTRACER_HEADER_ONLY=1 "];
        libArg =
          If[useLib,
            " '" <> libPath <> "'",
            ""];
        ntLog[
          "[time]   generator engine: ",
          If[useLib,
            "linking " <> libPath,
            "header-only (libNumTracer.a not found)"]];
        mainObj = bin <> "_main.o";
        unitObjs = Table[bin <> "_u" <> ToString[u - 1] <> ".o", {u, 1, Length[unitFiles]}];
(* Main-TU optimisation level. Compile and RUN each happen exactly ONCE per generation, so the
   only thing worth minimising is their SUM — an optimisation level is not "safer" for being
   higher, it is just a different split of the same total.

   -O2 IS DOMINATED and is no longer offered. Measured on the four-quark Fierz gate (non-dressed,
   100 diagrams, 16223 distinct traces), main TU only:
       -O0   compile  2.75 s   run 9.30 s   total 12.1 s
       -O1   compile 20.43 s   run 0.60 s   total 21.0 s
       -O2   compile 47.40 s   run 0.55 s   total 48.0 s
   so -O2 buys 0.05 s of run for 27 s of compile. The dressed path was already measured the same
   way (28.7 s vs 17.5 s compile, IDENTICAL run) — which is why the dressed default was -O1. The
   non-dressed default was -O2 only because that TU was assumed cheap to compile; that stopped
   being true once the per-net data tables dominated it, and it is what made the Fierz gate look
   like it hung (>600 s at -O2 before those tables were hash-consed).

   -O0 is NOT auto-selected. It wins on the Fierz gate (12.1 s vs 21.0 s), so an nSub-gated
   "-O0 for small flows" rule looks obviously right — and it is WRONG twice over.

   First, nSub does not predict the run. Tried with a cutoff of 20000 it made the codegen suite
   regenerate ZA4_147 for >18 minutes (from seconds), and that flow has nSub = 294 — twenty
   times FEWER distinct traces than Fierz's 3831, but each vastly bigger (dense four-gluon
   vertex over the full tensor basis). The run scales with trace SIZE, not count.

   Second, the real discriminator is DRESSED vs not: the DPoly contraction loop must stay
   inlined, the MPoly one is far less sensitive. Measured on ZAqbq (all_tensors, dressed,
   nSub = 306), main TU only, identical 2007 KB kernel from all three:
       -O0   compile 2.37 s   run 34.26 s   total 36.63 s
       -O1   compile 5.11 s   run  1.31 s   total  6.42 s
       -O2   compile 8.94 s   run  1.27 s   total 10.21 s
   (a third independent confirmation that -O2 is dominated).

   But "-O0 when not dressed" is fragile too: Fierz is non-dressed and STILL loses 16x of run
   at -O0 (0.60 s -> 9.92 s); it only wins because its compile saving happens to be larger, and
   a denser non-dressed flow would flip that. The asymmetry decides it — -O1 costs at most a
   bounded ~18 s of compile, while -O0 can cost 5.7x (ZAqbq) or minutes (ZA4_147), and a slow
   run looks exactly like a hang. So: -O1 always, with NT_GEN_MAIN_OPT=-O0 available for flows
   known to be small and non-dressed, where it is a genuine ~2x win. *)
        mainOpt =
          With[{e = Environment["NT_GEN_MAIN_OPT"]},
            Which[
              StringQ[e] && e =!= "",
                e,
              True,
                "-O1"]];
        ntLog["[time]   generator main TU: ", mainOpt, " (nSub = ", $ntGenNSub, "; NT_GEN_MAIN_OPT=-O0 is a large win on SMALL flows, but see the note above)"];
(* RAM-bounded parallel compile: run at most $ntCompileJobs `cxx` at once (xargs -P), and cap
   each at ~17 GB virtual (ulimit -v). Peak RAM ~ jobs x per-unit; with the unit count scaled so
   each TU is small (~12 net-builders) this stays well under the machine limit for any flow size.
   Both phases redirect compiler stdout+stderr to `clog` (compile truncates, link appends) so the
   genfail message can quote the actual g++ diagnostic, not just the exit code.

   -fno-exceptions -fno-rtti (UNIT TUs only): each net-builder unit is thousands of braced-init-
   lists of destructible temporaries (DiracNet / DChainTok / DSlot / NetVal); at -O0 emitting the
   exception-cleanup landing pads for them is ~90% of the compile and superlinear in unit size
   (a 188 KB unit measured 15.3 s -> 1.3 s, RSS 660 MB -> 280 MB, with exceptions off). The unit
   code never throws/catches, and the library's internal guards route through NT_THROW
   (core/config.hpp), which degrades to abort() there — so this only removes dead cleanup code.
   The MAIN TU keeps exceptions: it emits a try/catch thread-pool fallback (see the parWork
   template below, `catch(const std::system_error&)`), which is ill-formed under -fno-exceptions. *)
        pcmd = "printf '%s\\0' " <> StringRiffle[("\"" <> # <> "\"")& /@ Join[{"(ulimit -v 17000000; " <> cxx <> " -std=c++20 -ftemplate-depth=4000 " <> mainOpt <> hoDef <> "-pthread -I '" <> incDir <> "' -c '" <> genFile <> "' -o '" <> mainObj <> "')"}, Table["(ulimit -v 17000000; " <> cxx <> " -std=c++20 -ftemplate-depth=4000 -O0 -fno-exceptions -fno-rtti" <> hoDef <> "-I '" <> incDir <> "' -c '" <> unitFiles[[u]] <> "' -o '" <> unitObjs[[u]] <> "')", {u, 1, Length[unitFiles]}]], " "] <> " | xargs -0 -P " <> ToString[$ntCompileJobs] <> " -I CMD bash -c CMD > '" <> clog <> "' 2>&1";
        lcmd = cxx <> " -pthread '" <> mainObj <> "' " <> StringRiffle[("'" <> # <> "'")& /@ unitObjs, " "] <> libArg <> " -o '" <> bin <> "' >> '" <> clog <> "' 2>&1";
        {tcc, cc} =
          AbsoluteTiming[
            Run[pcmd];
            Run[lcmd]];
        Print["[time]   generator compile (", cxx, ", ", Length[unitFiles], " parallel units + main): ", tcc, " s"];
        If[cc =!= 0,
          Module[{
            lines =
              If[FileExistsQ[clog],
                StringSplit[ReadString[clog], "\n"],
                {}],
            head,
            nshow = 50},
            head = Take[lines, UpTo[nshow]];
            Message[
              mkGenerateKernel::genfail,
              cxx <> " compile/link rc=" <> ToString[cc] <> "\n--- first " <> ToString[Length[head]] <> " of " <> ToString[Length[lines]] <> " line(s) of " <> clog <> " ---\n" <> StringRiffle[head, "\n"] <>
                If[Length[lines] > nshow,
                  "\n... (" <> ToString[Length[lines] - nshow] <> " more line(s) in " <> clog <> ")",
                  ""]];
            Abort[]]]];
(* Free the per-generation codegen memo caches BEFORE launching the generator subprocess (hygiene:
   they're needed only to EMIT the source, already on disk, and re-cleared next generation anyway).
   NOTE (measured 2026-07-22): this is MINOR — the four caches together are only ~50 MB. The real
   WolframKernel resident tax that co-resides with the generator (~2 GB at 655 diagrams, ~3.7 GB at
   3350) is the `ntk` ITSELF (the front-end Diagrams/Components analysis, a function argument), which
   is not freeable here. The effective RAM lever is instead FEWER diagrams — propagator collection
   folds 3350 -> 655 and the Wolfram tax with it. *)
      If[$NumTracerVerbose, ntLog["[prof] pre-run  MemoryInUse=", Round[MemoryInUse[]/1048576.], " MB  RSS=", Round[ntWolframRssMB[]], " MB"]];
      $ctCache = <||>; $odCache = <||>; $dsCache = <||>; $dslCache = <||>;
      If[$NumTracerVerbose, ntLog["[prof] post-free MemoryInUse=", Round[MemoryInUse[]/1048576.], " MB  RSS=", Round[ntWolframRssMB[]], " MB"]];
(* run into a TEMP file, validate (rc==0 AND non-empty), then move into place — so a crashed
   generator (e.g. thread-limited Run[]) never silently truncates the committed header. *)
      Module[{tmp = headerFile <> ".tmp", rc, sz, trun},
        {trun, rc} =
          AbsoluteTiming[
            Run[
              "'" <> bin <> "' -n '" <> ns <> "' -d '" <> decor <> "'" <>
                If[OptionValue["FullParallel"],
                  " -p",
                  ""
                ] <> " > '" <> tmp <> "'"]];
        Print["[time]   generator run (reduce+rebase+lower): ", trun, " s"];
        sz =
          If[FileExistsQ[tmp],
            FileByteCount[tmp],
            0];
        If[rc =!= 0 || sz < 64,
          If[FileExistsQ[tmp],
            DeleteFile[tmp]];
          Message[mkGenerateKernel::genfail, "generator run rc=" <> ToString[rc] <> " bytes=" <> ToString[sz] <> " (committed header left intact)"];
          Abort[]];
        CopyFile[tmp, headerFile, OverwriteTarget -> True];
        DeleteFile[tmp];
        Print["wrote header: ", headerFile, " (", sz, " bytes)"]]];
(* semantic complexQ: with the real traces now generated, probe whether Im(integrand) actually
   vanishes (the projector-i / colour-i factors usually cancel to a real flow). If so, re-emit a
   REAL (double) kernel — losslessly, since Im≡0 — so the real DiFfRG integrators bind directly and
   no std::complex survives into device code. Only meaningful once the traces exist (RunGenerator). *)
    If[complexQ && TrueQ[OptionValue["RunGenerator"]] && TrueQ[OptionValue["RealProbe"]],
      Module[{
        verdict =
          numericImagProbeRealQ[
            integrand,
            args,
            fillArgs,
            angleDefs,
            angleDecls,
            nsHome,
            headerFile,
            $drTable,
            "TraceArrayDecl" ->
              If[crossCSE,
                tarrDecl,
                ""]],
        extraPriv = {}},
        Switch[verdict,
          "Pure", (* real value: project to Re. Wrap trace tokens in ntRe FIRST so the kernel is
provably double-typed even if a trace function is complex-typed — the probe
proved |full-proj|≈0, i.e. Σ Im(c)·tr ≈ 0, so dropping the ntIm terms is exact
and ntRe(double)=passthrough / ntRe(complex)=.real() keeps it real arithmetic. *)ntLog["[probe] real & projection exact -> clean REAL (double) kernel"];
            integrand = (integrand /. s_String :> Global`ntRe[s]) /. Complex[a_, b_] :> a;
            complexQ = False;
            extraPriv = {ntReImDefs[decor]},
          "RePart", (* real value but a complex trace: re/im split -> Σ[Re(c)·tr.real() − Im(c)·tr.imag()] *)ntLog["[probe] real value via complex trace(s) -> double kernel (re/im split)"];
(* The integrand is LINEAR in the trace tokens (strings). Re(Σ c·tr) = Σ[Re(c)·ntRe(tr) −
   Im(c)·ntIm(tr)]. The naive `(… /. s:>ntRe[s]+I ntIm[s]) /. Complex[a_,b_]:>a` is WRONG:
   Mathematica keeps `i·X·(ntRe+I ntIm)` as an UNEXPANDED product, so `/. Complex:>a` zeroes
   the leading `i` factor and DROPS the whole term — silently losing the real −X·ntIm(tr)
   contribution of every complex trace. Instead split each token's coefficient into real /
   imaginary parts via an `ii`-substitution (I → real symbol ii), which keeps the dressing
   coefficients FACTORED (unlike ComplexExpand, which un-factors them and defeats COEN's CSE). *)
            integrand =
              Module[{toks = Union[Cases[integrand, _String, Infinity]], tsym, lin, ii},
                tsym = AssociationThread[toks -> Table[Unique["tr$"], {Length[toks]}]];
                lin = integrand /. tsym;
                Total[
                    Function[t,
                        Module[{cc = Coefficient[lin, tsym[t]] /. Complex[ar_, ai_] :> ar + ii * ai},
                          Coefficient[cc, ii, 0] * Global`ntRe[t] - Coefficient[cc, ii, 1] * Global`ntIm[t]
                        ]
                      ] /@ toks
                  ] + ((lin /. Thread[Values[tsym] -> 0]) /. Complex[ar_, ai_] :> ar)];
            complexQ = False;
            extraPriv = {ntReImDefs[decor]},
          _,
            ntLog["[probe] imaginary part survives -> keeping the complex kernel (consumer takes Re)"]];
        If[verdict === "Pure" || verdict === "RePart",
          kernelFn = FunKit`MakeCppFunction[integrand, "Name" -> "kernel", "Prefix" -> decor, "Return" -> "auto", "CodeParser" -> "Cpp", "Parameters" -> kernelParams, "Body" -> preamble];
          classStr = ntKernelClass[name, {kernelFn, constFn}, decor, regTemplate, regAlias, extraPriv];
          header = FunKit`MakeCppHeader["Includes" -> Join[extraInc, ntRuntimeIncludes[runInc], {"numtracer/sun/sun_data.hpp", hdrInc}], "Body" -> ntWrapBody[kns, classStr, name]]
        ]]];
    (* kernel header (write-if-changed). *)
    If[FileExistsQ[kernelFile] && Import[kernelFile, "Text"] === header,
      Print["unchanged: ", kernelFile],
      ntExportCpp[kernelFile, header];
      Print["wrote kernel: ", kernelFile]];
    kernelFile];

(* ---- MakeNTKernel: the public kernel emitter. --------------------------------------
   MakeNTKernel[ntk, genFile, kernelFile, tracesFile] emits the numeric matrix-product kernel:
   a build-time generator program (genFile, + net-builder units + decl header), run to produce the
   committed straight-line traces header (tracesFile), and the kernel header (kernelFile) that fills
   the fundamental symbols and calls the traces. Options are forwarded to the generator
   (see Options[mkGenerateKernel] for the set). *)

Options[MakeNTKernel] = {"Name" -> "nt_kernel", "Namespace" -> Automatic, "Dressings" -> {}, "ScalarParams" -> {}, "ADParams" -> {}, "Decorator" -> "static inline", "IncludeDir" -> Automatic, "RunGenerator" -> True, "FullParallel" -> False, "AngleDefs" -> {}, "CrossTraceCSE" -> False, "GlobalCollect" -> True, "NumericContract" -> False, "Components" -> Automatic, "SymbolDefs" -> <||>, "RuntimeInclude" -> "numtracer/codegen/runtime.hpp", "ExtraIncludes" -> {}, "KernelNamespace" -> "numtracer_kernels", "SupportNamespace" -> "numtracer", "DressingType" -> Automatic, "RegulatorTemplate" -> False, "RegulatorAlias" -> False, "RealProbe" -> True, "PruneRealTraces" -> False, "Constant" -> 0.};

MakeNTKernel::disconnectmix = "Diagram `1` disconnects into >= 2 Dirac/colour trace components (a product of independent Dirac traces, a genuine >=2-loop structure). The numeric backend handles a single Dirac/colour trace times any number of disconnected pure-Lorentz scalars (factored), but does not yet multiply two or more independent Dirac traces.";

MakeNTKernel::nfiles = "MakeNTKernel needs three output files: MakeNTKernel[ntk, genFile, kernelFile, tracesFile].";

(* Explicit opts win (first match); MakeNTKernel's own defaults carry through anything not passed —
   so `SetOptions[MakeNTKernel, ...]` (used by the test setups to opt into the shim + the DiFfRG
   namespace) propagates to the generator without editing every call site. *)

MakeNTKernel[ntk : NTKernel[_], file_, opts : OptionsPattern[]] := (
    Message[MakeNTKernel::nfiles];
    Abort[]);

MakeNTKernel[ntk : NTKernel[_], genFile_, kernelFile_, tracesFile_, opts : OptionsPattern[]] :=
  mkGenerateKernel[ntk, genFile, kernelFile, tracesFile, "NumericContract" -> True, Sequence @@ FilterRules[Join[{opts}, Options[MakeNTKernel]], Options[mkGenerateKernel]]];
