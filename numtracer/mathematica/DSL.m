(* ::Package:: *)

(* DSL analysis: classify the heads, split each term into independent contraction
   components, and allocate the env-id layout the `et` builders consume.

   Loaded inside NumTracer`Private` by NumTracer.m — public symbols (NumTrace, the
   nt* heads) already exist in the NumTracer` context. *)

(* ---- head classification ---------------------------------------------------- *)

(* A factor that participates in the tensor contraction (vs. a scalar coefficient).
   ntEpsilon (the Levi-Civita ε_{μνρσ}) is a Lorentz tensor produced by the γ5 trace theorem
   (gammaTraceSum5 in Codegen.m); it carries four Lorentz labels like a metric pair. *)
(* ntVec[q, mu] with a SYMBOLIC label is a tensor leg; ntVec[q, i_Integer] is the scalar component
   q_i (0-based, 0=temporal), resolved by the frame in the coeff like ntSP — so it is NOT a tensor.
   Express that by matching only a non-integer second arg here (an integer arg falls through to False). *)
tensorQ[_ntMetric | ntVec[_, Except[_Integer]] | _ntTransProj | _ntLongProj |
        _ntElectricProj | _ntMagneticProj | _ntSUNf | _ntSUNDeltaAdj |
        _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac | _ntSUNT | _ntSUNDeltaFund | _ntEpsilon |
        _ntSUNDiagFund | _ntSUNDiagAdj | _ntDressedNum] = True;
tensorQ[_] = False;

(* The 4 SU(N) group heads carry their rank N as the FIRST argument, so colour SU(Nc),
   flavour/isospin SU(Nf), and any further group coexist in one network at the right rank.
   Lorentz, Dirac, and the SU(N) groups contract in disjoint index spaces (matched by shared
   label), so the label sets are treated uniformly while the builders stay distinct. *)
(* ntSUNDiag{Fund,Adj}[N, i, j, spec, scale] — a group δ that carries a PER-COMPONENT dressing.
   `spec` is a rules list {c -> name, ..., Default -> defName} of 1-based component indices to
   distinctly-named scalar dressing symbols (evaluated at the kinematic `scale`); unnamed components
   collapse to Default when given, else drop. They classify exactly like the matching plain δ
   (fundamental / adjoint) for index bookkeeping; the per-component dressing is folded numerically
   by sun_value_dressed at codegen time. *)
adjointSUNQ[_ntSUNf | _ntSUNDeltaAdj | _ntSUNDiagAdj] = True;  (* group-adjoint heads (bridge with Lorentz) *)
adjointSUNQ[_] = False;
fundamentalSUNQ[_ntSUNT | _ntSUNDeltaFund | _ntSUNDiagFund] = True;   (* group-fundamental heads (quark-line) *)
fundamentalSUNQ[_] = False;

(* The SU(N) rank a group head builds against: its leading argument. *)
sunRankOf[h_] := First[h];

(* The contraction index labels carried by a tensor factor (momentum arg dropped). *)
labelsOf[ntMetric[mu_, nu_]]         := {mu, nu};
labelsOf[ntVec[_, mu_]]              := {mu};
labelsOf[ntTransProj[_, mu_, nu_]]   := {mu, nu};
labelsOf[ntLongProj[_, mu_, nu_]]    := {mu, nu};
labelsOf[ntElectricProj[_, mu_, nu_]] := {mu, nu};
labelsOf[ntMagneticProj[_, mu_, nu_]] := {mu, nu};
labelsOf[ntSUNf[_, a_, b_, c_]]      := {a, b, c};
labelsOf[ntSUNDeltaAdj[_, a_, b_]]   := {a, b};
labelsOf[ntGamma[mu_, din_, dout_]]  := {mu, din, dout};
labelsOf[ntGamma5[din_, dout_]]      := {din, dout};
(* ntSigma's OPEN Lorentz labels are its FREE legs (gluon ids); slashed legs carry a momentum, not an
   open id. The spinor axes din,dout are always open. *)
labelsOf[ntSigma[legA_, legB_, din_, dout_]] :=
  Join[Cases[{legA, legB}, {"free", mu_} :> mu], {din, dout}];
labelsOf[ntDeltaDirac[din_, dout_]]  := {din, dout};
labelsOf[ntSUNT[_, a_, i_, j_]]      := {a, i, j};
labelsOf[ntSUNDeltaFund[_, i_, j_]]  := {i, j};
labelsOf[ntSUNDiagFund[_, i_, j_, _, _]] := {i, j};  (* per-component-dressed fundamental δ *)
labelsOf[ntSUNDiagAdj[_, a_, b_, _, _]]  := {a, b};  (* per-component-dressed adjoint δ *)
labelsOf[ntEpsilon[a_, b_, c_, d_]]  := {a, b, c, d};
(* ntDressedNum[options, din, dout] — a dressed propagator numerator kept EAGER (symbolic dressing
   collection). `options` is a tensor-FREE list of {coeffExpr, spec} (spec = {"ident"} | {"slash",vlc})
   so it exposes only its spinor in/out labels; the dressing collection is folded in C++ (one DPoly
   trace) instead of distributing the numerator into 2^D diagrams. *)
labelsOf[ntDressedNum[_, din_, dout_]] := {din, dout};

(* Spinor (Dirac) axis labels of a head — the subset of labelsOf that lives in the
   spinor index space. Used to give spinor axes a disjoint id range so the et engine
   never contracts a spinor axis against a Lorentz/colour axis sharing an id. *)
spinorLabelsHead[ntGamma[_, din_, dout_]] := {din, dout};
spinorLabelsHead[ntGamma5[din_, dout_]]   := {din, dout};
spinorLabelsHead[ntSigma[_, _, din_, dout_]] := {din, dout};
spinorLabelsHead[ntDeltaDirac[din_, dout_]] := {din, dout};
spinorLabelsHead[ntDressedNum[_, din_, dout_]] := {din, dout};
spinorLabelsHead[_]                       := {};
(* All spinor labels anywhere under a (sub)expression. *)
allSpinorLabels[e_] := DeleteDuplicates @ Flatten @ Cases[e, h_?tensorQ :> spinorLabelsHead[h], {0, Infinity}];

(* The momentum a factor carries (the source of runtime Var leaves), or None. *)
momentumOf[ntVec[q_, _]]             := q;
momentumOf[ntTransProj[q_, _, _]]    := q;
momentumOf[ntLongProj[q_, _, _]]     := q;
momentumOf[ntElectricProj[q_, _, _]] := q;
momentumOf[ntMagneticProj[q_, _, _]] := q;
momentumOf[_]                        := None;

(* Whether a momentum needs a 1/q^2 env slot (it sits inside a projector with a full denominator). *)
needsInvQ[ntTransProj[__]]    = True;
needsInvQ[ntLongProj[__]]     = True;
needsInvQ[ntElectricProj[__]] = True;  (* electric P_E = P_T − P_M uses 1/q² (and 1/|q⃗|², see needsInvSQ) *)
needsInvQ[_]                  = False;

(* Whether a momentum needs a SPATIAL 1/|q⃗|² env slot (finite-T electric/magnetic projectors). *)
needsInvSQ[ntElectricProj[__]] = True;
needsInvSQ[ntMagneticProj[__]] = True;
needsInvSQ[_]                  = False;

(* ---- self-trace normalization ----------------------------------------------- *)

(* A tensor with a repeated index is a self-trace (e.g. P^mu_mu). The et engine
   contracts pairwise BETWEEN tensors and never self-contracts one, so we relabel
   the second occurrence and insert the matching identity (Lorentz metric for a
   Lorentz index, adjoint delta for a colour index): P^mu_mu = P^{mu nu} d_{mu nu}.
   This keeps every index appearing on two distinct tensors, as the engine needs. *)
splitSelfTraces[factors_List] := Module[{res = {}, conns = {}},
  Function[f, If[! tensorQ[f], AppendTo[res, f],
    (* The connecting identity reuses the SAME group rank N as the head it closes: an
       adjoint group self-trace closes with ntSUNDeltaAdj[N,..], a fundamental one with
       ntSUNDeltaFund[N,..], a Lorentz/Dirac one with the metric. *)
    Module[{dups, relabeled = f, conn = Which[
        adjointSUNQ[f],     With[{n = sunRankOf[f]}, ntSUNDeltaAdj[n, ##] &],
        fundamentalSUNQ[f], With[{n = sunRankOf[f]}, ntSUNDeltaFund[n, ##] &],
        True,               ntMetric]},                          (* Lorentz/Dirac (unchanged) *)
      dups = Cases[Tally[labelsOf[f]], {l_, c_} /; c >= 2 :> l];
      Do[With[{fresh = Unique["st"]},
           relabeled = ReplacePart[relabeled, Last[Position[relabeled, l, {1}]] -> fresh];  (* relabel 2nd occurrence *)
           AppendTo[conns, conn[l, fresh]]], {l, dups}];
      AppendTo[res, relabeled]]]] /@ factors;
  Join[res, conns]
];

(* ---- free indices & scalar test (work through Plus/Times, for eager summation) ---- *)

(* Whether a (sub)expression carries no tensor head — a pure scalar coefficient. *)
scalarQ[e_] := FreeQ[e, _ntMetric | ntVec[_, Except[_Integer]] | _ntTransProj | _ntLongProj |
                       _ntElectricProj | _ntMagneticProj | _ntSUNf | _ntSUNDeltaAdj |
                       _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac | _ntSUNT | _ntSUNDeltaFund | _ntEpsilon |
                       _ntSUNDiagFund | _ntSUNDiagAdj | _ntDressedNum];

(* The free (uncontracted) index labels of a tensor (sub)expression. A product sums
   indices that appear twice (free = appear once); a sum's summands share free indices
   (a vertex's legs); a bare head exposes all its labels. Used to group components and
   to align et::add operands. *)
freeIdx[e_] := Which[
  tensorQ[e],        labelsOf[e],
  Head[e] === Plus,  freeIdx[First[List @@ e]],
  (* A tensor^n is n copies sharing the SAME labels — a closed self-contraction (see compileTInv),
     so it exposes NO free index. Spelled out rather than left to the True branch below, which
     would return {} for the wrong reason and hide a malformed Power. *)
  Head[e] === Power && IntegerQ[e[[2]]] && e[[2]] >= 2 && ! scalarQ[e], {},
  (* count == 1 is free; count == 2 is contracted. checkLabels guarantees no label occurs more
     often, so the old OddQ test was equivalent — but OddQ would silently call a 4x label
     "contracted" and a 3x label "free" if the guard were ever bypassed. *)
  Head[e] === Times, Cases[Tally[Flatten[freeIdx /@ (List @@ e)]], {l_, c_} /; c == 1 :> l],
  True,              {}];

(* All index labels anywhere under a (sub)expression (free + internally summed) — the
   set that needs an axis id. *)
allLabels[e_] := DeleteDuplicates @ Flatten @ Cases[e, h_?tensorQ :> labelsOf[h], {0, Infinity}];

(* ---- connected components over shared free indices -------------------------- *)

(* Group tensor factors (heads OR Plus-vertices) into independent sub-networks: two
   factors are connected if their free indices intersect. Each component contracts to
   a scalar (colour stays a separate constant); keeps every contraction small. *)
connectedComponents[factors_List] := Module[{verts, edges, graph},
  verts = Range[Length[factors]];
  edges = Select[Subsets[verts, {2}],
    IntersectingQ[freeIdx[factors[[#[[1]]]]], freeIdx[factors[[#[[2]]]]]] &];
  graph = Graph[verts, UndirectedEdge @@@ edges];
  factors[[#]] & /@ (ConnectedComponents[graph])
];

(* Greedy contraction order: keep each successive factor sharing a free index with the
   running set, so contract_all's intermediates stay low-rank (never outer-product the
   whole thing). *)
orderFactors[fs_List] := Module[{rem = fs, out, used},
  If[fs === {}, Return[{}]];
  out = {First[rem]}; used = freeIdx[First[rem]]; rem = Rest[rem];
  While[rem =!= {},
    With[{pick = SelectFirst[rem, IntersectingQ[freeIdx[#], used] &, First[rem]]},
      AppendTo[out, pick]; used = Union[used, freeIdx[pick]];
      rem = DeleteCases[rem, pick, {1}, 1]]];
  out];

(* ---- sector-bridge expansion (keep colour and Lorentz contractions separate) ---- *)

(* Eager summation (et::add) is a win only when a structure-sum lives in ONE sector — a sum of
   Lorentz structures with a single colour factor pulled out (e.g. a 3-gluon vertex) stays one
   small Lorentz tensor. A SECTOR-BRIDGING sum is different: its colour sum is *correlated* with
   its Lorentz sum (each term pairs one colour structure with one Lorentz structure, e.g. a
   4-gluon vertex's f.f ⊗ metric.metric). et::add-ing those fuses the colour (dim N^2-1) and
   Lorentz (dim 4) axes into one ETensor whose entry count (and type) explodes — building it alone
   times out. The cure: DISTRIBUTE such a sum, so each term is colour/Lorentz index-disjoint again
   — colour folds to a constant, Lorentz stays small. That is linear in the (few) structures, not
   the product-of-sums blow-up. *)
sectorBridgeQ[p_] := (! FreeQ[p, _ntSUNf | _ntSUNDeltaAdj | _ntSUNDiagAdj]) &&
                     (! FreeQ[p, _ntMetric | ntVec[_, Except[_Integer]] | _ntTransProj | _ntLongProj |
                                 _ntElectricProj | _ntMagneticProj]);

(* The scalar (non-tensor) coefficient of a single summand. *)
scalarCoeffOf[t_] := Times @@ Select[If[Head[t] === Times, List @@ t, {t}], scalarQ];
(* A TENSOR-structure sum whose per-summand scalar coefficients are not all numeric must be
   DISTRIBUTED too: the et engine's eager add (et::add) can only fold NUMERIC per-structure
   scalars (it scales each structure by a compile-time Lit). A fermion propagator numerator
   Mq*deltaDirac + (.. dressings ..)*(gamma.vec) carries runtime dressing coefficients, so
   each Dirac structure must become its own diagram with that dressing as a scalar coeff.
   (Pure-scalar dressing sums — no tensor head — are left intact as coefficients.) *)
dressedStructureSumQ[p_Plus] := (! scalarQ[p]) && AnyTrue[List @@ p, ! NumericQ[scalarCoeffOf[#]] &];
(* Symbolic dressing collection: when $ntDressCollect is True a
   DIRAC dressed structure sum (a propagator numerator like Mq·δ + Z(p)·γ·p, all terms sharing the
   spinor in/out indices) is NOT distributed — it is kept eager and rewritten to a single ntDressedNum
   slot (rewriteDressedNums), so the diagram folds to ONE DPoly trace carrying its dressings as
   `dress` env leaves instead of exploding into 2^D diagrams. A colour↔Lorentz sectorBridge sum still
   distributes (sectors must split), and a non-Dirac dressed sum (no shared spinor indices) still
   distributes (the collection only handles the Dirac numerator case here). Default: $ntDressCollect
   False = today's full distribution. *)
$ntDressCollect = False;
(* The OPEN (free) spinor labels of a summand: spinor indices appearing an odd number of times. *)
openSpinorOf[t_] := Cases[Tally[Flatten[Cases[t, h_?tensorQ :> spinorLabelsHead[h], {0, Infinity}]]],
  {l_, c_} /; OddQ[c] :> l];
(* A propagator-numerator structure sum: every term carries Dirac structure and the SAME pair of open
   spinor indices (so it connects one spinor-in to one spinor-out, like Mq·δ[a,b] + Z·γ[μ,a,b]·vec[μ]). *)
diracNumeratorSumQ[p_Plus] := Module[{terms = List @@ p, opens},
  opens = openSpinorOf /@ terms;
  AllTrue[terms, ! FreeQ[#, _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac] &] &&
    AllTrue[opens, Length[#] === 2 &] && SameQ @@ (Sort /@ opens)];
(* A dressed Dirac sum is collectible only if every term is a propagator numerator we can ACTUALLY
   decompose (ident / slash, no open gluon leg). diracNumeratorSumQ is the cheap pre-filter (same 2
   open spinor indices, Dirac structure); dressedNumDecompose =!= $Failed is the exact test — a sum
   with an open Lorentz leg (a multi-structure quark-gluon / σ vertex) fails it and stays on the
   distribution path (so it never survives un-distributed AND un-rewritten → no eagernn abort). *)
collectibleDiracSumQ[p_Plus] := dressedStructureSumQ[p] && ! sectorBridgeQ[p] &&
  diracNumeratorSumQ[p] && dressedNumDecompose[p] =!= $Failed;
collectibleDiracSumQ[_] := False;
distributeQ[p_] := sectorBridgeQ[p] ||
  (dressedStructureSumQ[p] && ! (TrueQ[$ntDressCollect] && collectibleDiracSumQ[p]));

(* The SET of γ-count parities the diagram's terms would carry IF every EAGER Dirac Plus were
   distributed — {0}, {1}, or {0,1}. Used for the odd-trace-vanishing verdict: a closed trace is
   identically zero only when EVERY branch has an odd γ count, i.e. the set is exactly {1}.
   A diagram-global `Count[_ntGamma, Infinity]` cannot decide this — it sums γ's ACROSS the branches
   of an eager Plus and so reports a parity no actual term has (e.g. a projector 1 + p̸ against two
   γ's counts 3 → "odd" → the whole diagram is dropped, taking its non-vanishing even branch with it).
   Computed WITHOUT Expand: the sets are capped at two elements and fold pairwise, so this is linear
   in the factor count rather than a product-of-sums blow-up.
   Only _ntGamma is counted, matching the old test: ntSigma (2 γ) and ntDeltaDirac (0 γ) are both
   EVEN and so correctly contribute parity 0. ntDressedNum cannot appear yet (rewriteDressedNums runs
   later, inside analyseDiagram). *)
diracParities[e_Plus] := Union @@ (diracParities /@ (List @@ e));
diracParities[e_Times] := Fold[Union[Flatten[Mod[Outer[Plus, #1, #2], 2]]] &,
                               {0}, diracParities /@ (List @@ e)];
(* b^n is n copies of b multiplied. Do NOT write this as diracParities[Times @@ ConstantArray[b, n]]:
   Times immediately re-collapses the copies back to b^n, so that recurses until $RecursionLimit and
   returns garbage — and it fires on ORDINARY SCALAR powers (sp[p2,p3]^2 in a Gram determinant), not
   just tensor powers. Fold the parity arithmetically instead: n copies of a definite parity p give
   n*p mod 2, and if b itself has both parities then so does any power of it. *)
diracParities[Power[b_, n_Integer?Positive]] := With[{s = diracParities[b]},
  If[Length[s] > 1, {0, 1}, {Mod[n * First[s], 2]}]];
(* {0, Infinity}, NOT Infinity: this catch-all is reached with a BARE factor (a single ntGamma head
   pulled out of the enclosing Times), and level spec Infinity means {1, Infinity} — it would skip
   level 0 and count that γ as zero. The old diagram-global test never hit this because it always
   ran on a whole Times, where every γ sits at level 1. *)
diracParities[e_] := {Mod[Count[e, _ntGamma, {0, Infinity}], 2]};

(* The odd-trace verdict itself: a γ5-free diagram every one of whose branches is an odd closed
   γ trace is identically zero and can be pruned before any trace runs. *)
vanishingOddTraceQ[diagram_] := FreeQ[diagram, _ntGamma5] && diracParities[diagram] === {1};

(* The slash momentum of a γ·p term: the {coeff, q} pairs of the ntVec factors sharing the γ's Lorentz
   index μ (a signed linear combination like p−l, exactly like ntSigmaLeg). $Failed if ill-formed. *)
dressedVlc[facs_, mu_] := Module[{lf = Times @@ Select[facs, ! FreeQ[#, ntVec[_, mu]] &], terms},
  If[lf === 1, Return[$Failed]];
  terms = If[Head[lf] === Plus, List @@ lf, {lf}];
  Catch[Function[trm, Module[{q = Cases[trm, ntVec[qq_, m_] /; m === mu :> qq, {0, Infinity}],
                              c = trm /. ntVec[_, m_] /; m === mu :> 1},
     If[Length[q] =!= 1 || ! FreeQ[c, ntVec] || ! NumericQ[c], Throw[$Failed, "dvlc"]];
     {c, First[q]}]] /@ terms, "dvlc"]];

(* Decompose one term of a dressed numerator into {din, dout, scalar, spec, otherTensors}: the Dirac
   in/out spinor labels; the scalar (dressing × numeric × denominator) coefficient; the Dirac structure
   spec — {"ident"} (spinor-δ) or {"slash", {{c,q}…}} (γ·p̸); and the NON-Dirac tensor factors (a colour
   δ on the quark line, etc.) that must be COMMON across the sum's terms and factor out. $Failed for an
   unsupported structure (σ, γ5, >1 γ). *)
dressedNumTerm[t_] := Module[{facs = If[Head[t] === Times, List @@ t, {t}], dt, scal, mu, vlc, other},
  scal = Times @@ Select[facs, scalarQ];
  dt = Select[facs, MatchQ[#, _ntDeltaDirac | _ntGamma] &];
  Which[
    MatchQ[dt, {ntDeltaDirac[_, _]}],
      other = Select[DeleteCases[facs, dt[[1]]], tensorQ];
      {dt[[1, 1]], dt[[1, 2]], scal, {"ident"}, Sort[other]},
    MatchQ[dt, {ntGamma[_, _, _]}],
      mu = dt[[1, 1]]; vlc = dressedVlc[facs, mu];
      If[vlc === $Failed, $Failed,
        (other = Select[DeleteCases[facs, dt[[1]]], tensorQ[#] && FreeQ[#, ntVec[_, mu]] &];
         {dt[[1, 2]], dt[[1, 3]], scal, {"slash", vlc}, Sort[other]})],
    True, $Failed]];

(* Multiset of multiplicative factors common to every list in `factLists` (min multiplicity). *)
commonFactorMultiset[factLists_] := Module[{cnts = Counts /@ factLists, keys},
  keys = Intersection @@ (Keys /@ cnts);
  Flatten[Function[k, ConstantArray[k, Min[(Lookup[#, k, 0] &) /@ cnts]]] /@ keys]];

dressedNumDecompose[p_Plus] := Module[
  {rows = dressedNumTerm /@ (List @@ p), din, dout, others, common, scalFacs, commonScal, opts},
  If[MemberQ[rows, $Failed], Return[$Failed]];
  {din, dout} = rows[[1, {1, 2}]];
  If[! AllTrue[rows, #[[1]] === din && #[[2]] === dout &], Return[$Failed]]; (* all terms din→dout *)
  others = rows[[All, 5]];
  If[! AllTrue[others, Sort[#] === Sort[others[[1]]] &], Return[$Failed]];   (* common colour factors *)
  common = others[[1]];
  (* factor the scalar coefficient common to every term (the propagator denominator, a flavour δ, …)
     OUT of the sum so it multiplies the whole numerator; only the per-structure residual (the
     dressing that differs: Mq vs Zq·…) stays inside the ntDressedNum options. *)
  scalFacs = Function[s, If[Head[s] === Times, List @@ s, {s}]][#[[3]]] & /@ rows;
  commonScal = commonFactorMultiset[scalFacs];
  opts = MapThread[Function[{r, sf},
     {Times @@ Fold[DeleteCases[#1, #2, {1}, 1] &, sf, commonScal], r[[4]]}], {rows, scalFacs}];
  (Times @@ common) * (Times @@ commonScal) * ntDressedNum[opts, din, dout]];

(* Rewrite each collectible dressed Dirac numerator (a surviving Plus factor) into one ntDressedNum
   token, and SPLICE the factored-out common colour / flavour / denominator product back into the
   factor list as separate factors — so the scalars (flavour δ, 1/denom) land in the diagram coeff
   (where contractFlavour can collapse a flavour chain) and the colour δ folds as its own factor. A
   factor that decomposes to $Failed is left as-is. Only active under $ntDressCollect. *)
rewriteDressedNums[factors_List] := If[! TrueQ[$ntDressCollect], factors,
  Flatten[Function[f, If[Head[f] === Plus && collectibleDiracSumQ[f],
     With[{r = dressedNumDecompose[f]},
       If[r === $Failed, {f}, If[Head[r] === Times, List @@ r, {r}]]], {f}]] /@ factors]];

(* Inverse of the ntDressedNum rewrite: expand a collected numerator back into its distributed Dirac
   structure sum (Σ coeff_i · {ntDeltaDirac | ntGamma·Σc·ntVec}). Used when re-distributing a
   collected diagram (the small-D gate in NumTrace). Each slash gets a fresh Lorentz label for its
   contracted leg. *)
expandDressedNum[ntDressedNum[opts_, din_, dout_]] := Plus @@ (Function[opt,
  opt[[1]] * Switch[opt[[2, 1]],
    "ident", ntDeltaDirac[din, dout],
    "slash", With[{mu = Unique["dexp"]},
       ntGamma[mu, din, dout] * (Plus @@ ((#[[1]] * ntVec[#[[2]], mu]) & /@ opt[[2, 2]]))]]] /@ opts);

(* Re-distribute ONE analysed (collected) diagram back to the non-collected path: rebuild its net,
   expand every ntDressedNum into the Dirac structure sum, distribute, drop odd-gamma traces, and
   re-analyse (with collection OFF so no ntDressedNum is recreated). Used by the small-D gate in
   NumTrace. *)
redistDiagram[diag_] := Block[{$ntDressCollect = False},
  Module[{net = diag["Coeff"] * Times @@ Flatten[(#["Factors"] &) /@ diag["Components"]]},
    net = net /. nd_ntDressedNum :> expandDressedNum[nd];
    With[{ex = expandBridges[net]},
      analyseDiagram /@ Select[If[Head[ex] === Plus, List @@ ex, {ex}],
        ! vanishingOddTraceQ[#] &]]]];

(* Distribute every colour<->Lorentz-bridging sum into its surrounding product (only that
   sum; single-sector sums are left intact for et::add). Turns a bridging diagram into a
   small linear sum of sector-separable diagrams. Explicit recursion (not a //. rule over
   Orderless Times, which backtracks catastrophically on a large net). *)
expandBridges[e_Plus] := Plus @@ (expandBridges /@ (List @@ e));
expandBridges[e_Times] := Module[{factors = List @@ e, bridge},
  bridge = FirstCase[factors, p_Plus /; distributeQ[p], Missing[]];
  If[MissingQ[bridge],
    Times @@ (expandBridges /@ factors),
    expandBridges[Plus @@ (Times @@ Append[DeleteCases[factors, bridge, {1}, 1], #] & /@ (List @@ bridge))]]];
(* Times is Flat+Orderless and auto-collects identical factors, so two byte-identical bridging sums
   become Power[sum, 2] BEFORE we ever see them. The selector above matches p_Plus only, so such a
   Power is never chosen as a bridge and would fall through the catch-all below UN-DISTRIBUTED —
   fusing the colour and Lorentz axes into one giant ETensor (the exact blow-up expandBridges
   exists to prevent), or worse being read by compileTInv as a closed self-contraction. It cannot
   simply be re-expanded here: Times would immediately re-collapse the copies unless their labels
   were freshened first, and two vertices that legitimately share every label are themselves a bug.
   No known flow produces this; fail loudly if one ever does. *)
NumTrace::bridgepow = "expandBridges: a colour<->Lorentz-bridging sum appears raised to the power \
`1`, i.e. as `1` byte-identical factors sharing every index label. Such a sum cannot be distributed \
(and identical labels on two distinct vertices are themselves malformed). Offending base:\n`2`";
expandBridges[e : Power[b_Plus, n_Integer]] /; n >= 2 && distributeQ[b] :=
  (Message[NumTrace::bridgepow, n, Short[b, 6]]; Abort[]);
expandBridges[e_] := e;

(* ---- env-id layout ---------------------------------------------------------- *)

(* Assign each distinct momentum a Base (4 consecutive Var ids) and, where a
   projector needs it, an Inv id holding 1/q^2. Var ids index the runtime renv[]. *)
buildEnv[momenta_List, invMomenta_List, invSMomenta_List] := Module[{env = <||>, base = 0, inv},
  Do[env[q] = <|"Base" -> base, "Inv" -> None, "InvS" -> None|>; base += 4, {q, momenta}];
  inv = base;
  Do[env[q]["Inv"] = inv++, {q, invMomenta}];      (* full 1/q² slots *)
  Do[env[q]["InvS"] = inv++, {q, invSMomenta}];    (* spatial 1/|q⃗|² slots (finite-T E/M projectors) *)
  {env, inv}  (* inv is now the total env size NEnv *)
];

(* Per-momentum component mask: bit i set <=> component i is structurally nonzero in
   the frame. Drives the et builder Mask template arg (zero components prune away). *)
frameMask[components_List] := FromDigits[Reverse[Boole[# =!= 0 && # =!= 0.] & /@ components], 2];

(* ---- NumTrace --------------------------------------------------------------- *)

(* Each SU(N) group head carries its own rank N as the first argument (baked in when the
   network is built — Global`Nc for colour, the FromFunKit "FlavourGroup" option for the
   isospin group), so NumTrace itself takes no group option. *)
Options[NumTrace] = {"Frame" -> <||>, "Args" -> {}, "Dressings" -> {}, "DressingCollection" -> True};

NumTrace[net_, OptionsPattern[]] := Module[
  {frame, args, dress, badRanks, diagrams, allMom, invMom, invSMom, env, nenv, diags},
  frame  = OptionValue["Frame"];
  args   = OptionValue["Args"];
  dress  = OptionValue["Dressings"];
  (* symbolic dressing collection: keep dressed Dirac numerators eager (one DPoly trace) instead of
     distributing into 2^D diagrams. Default False = today's behaviour. Set here so expandBridges
     (distributeQ) and analyseDiagram (rewriteDressedNums) both see it; each NumTrace call sets it
     from its option (default False), so it never leaks across kernels. *)
  $ntDressCollect = TrueQ[OptionValue["DressingCollection"]];

  (* SU(N) ranks must be compile-time integers (they pick the correct-dimension typed-out
     group matrices). Every group head's leading argument N is checked up front, so an
     undefined symbol aborts with a clear message rather than being baked into the kernel. *)
  badRanks = DeleteDuplicates @ Select[
    Cases[net, h : (_ntSUNf | _ntSUNDeltaAdj | _ntSUNT | _ntSUNDeltaFund | _ntSUNDiagFund | _ntSUNDiagAdj) :> sunRankOf[h], Infinity],
    ! (IntegerQ[#] && # >= 1) &];
  If[badRanks =!= {},
    Print["NumTrace: every SU(N) head must carry an integer rank N >= 1 as its first argument ",
          "(got ", badRanks, "). Set Global`Nc (SetNc[3]) before tracing colour heads, and pass ",
          "\"FlavourGroup\" -> n to FromFunKit for the isospin group."]; Abort[]];
  (* A closed quark flavour loop folds to Global`Nf. When SetNf[n] has been called Global`Nf is
     an integer and has already folded into the coefficients (nothing to check). Only when it is
     still the unset symbol AND it appears in the network is the flavour count genuinely undefined
     — check for the symbol then (guarding so the integer value is not mistaken for a literal). *)
  If[! IntegerQ[Global`Nf] && ! FreeQ[net, Global`Nf],
    Print["NumTrace: flavour Nf is not a defined integer (symbolic Nf present in the network). ",
          "Call SetNf[2] (TensorBases) before generating."]; Abort[]];

  (* the top-level sum is the (linear) sum of DIAGRAMS; keep single-sector vertex sums eager
     via et::add, but distribute colour<->Lorentz-bridging sums so the two sectors never fuse
     into one giant ETensor. *)
  ntLog["[prof] NumTrace expandBridges: ", First@AbsoluteTiming[
  diagrams = With[{ex = expandBridges[net]}, If[Head[ex] === Plus, List @@ ex, {ex}]];], " s"];
  (* odd-trace vanishing: a closed Dirac trace of an ODD number of gammas is identically zero
     (and would otherwise carry a spurious imaginary I^odd coefficient from the vertex/basis
     normalizations). Drop such diagrams. (gamma5-bearing traces are exempt from the rule.)
     The verdict is PER BRANCH (diracParities), not on a diagram-global gamma count: an EAGER Dirac
     Plus left standing by expandBridges can carry a VARIABLE gamma parity, and summing the count
     across its branches yields a parity no actual term has. Two ways that bites:
       * a COLLECTIBLE dressed propagator sum (Mq·δ + Z·γ·p̸, kept eager under collection) — its δ
         ("ident") branch has 0 gammas, its γ·p̸ ("slash") branch 1. (e.g. struct-4 `l̸1 γ^ρ` + a
         propagator counts 3 γ → "odd" → the diagram was dropped, killing its non-zero Mq·δ branch
         — the 1/4/7 collection bug.)
       * a MULTI-TERM PROJECTOR's Plus — all-numeric coefficients and single-sector, so distributeQ
         is False (left eager) AND collectibleDiracSumQ is False (that predicate needs a NON-numeric
         coefficient). It was exempted by neither, so the oblique-metric dual projector for e.g.
         AqbqDirect8 structure 7 traced to an identically-zero kernel with no error raised.
     Only an ALL-odd-branch diagram is dropped; a mixed-parity one is KEPT, and the per-branch
     matrix-product trace zeroes its odd-gamma branches on its own (splitColourGroupsInv expands any
     γ-bearing Plus and compiles each branch separately, so this filter is purely a pruning
     optimisation sitting upstream of already-correct code). *)
  diagrams = Select[diagrams, ! vanishingOddTraceQ[#] &];

  (* Validate every distributed diagram BEFORE analyseDiagram assigns axis ids: it gives ONE id
     per DISTINCT label (see below), so a label occurring 4x becomes four axes sharing an id and
     the et engine mis-pairs them into a silently wrong number. Checked per diagram (not on the
     raw net): only after expandBridges is each diagram a flat Times in which "1 = free,
     2 = contracted" is the actual invariant. *)
  ntLog["[prof] NumTrace checkLabels: ", First@AbsoluteTiming[
  With[{frees = MapIndexed[checkLabels[#1, First[#2]] &, diagrams]},
    ntLog["[labels] ", Length[diagrams], " diagram(s) validated; free-index set(s) = ",
      DeleteDuplicates[Sort /@ frees]]];], " s"];

  (* global env layout: every distinct momentum, and which ones need a 1/q^2 slot *)
  allMom = DeleteDuplicates @ Cases[net, f_?tensorQ :> momentumOf[f], Infinity] // DeleteCases[None];
  invMom = DeleteDuplicates @ Cases[net, f_?(needsInvQ) :> momentumOf[f], Infinity];
  invSMom = DeleteDuplicates @ Cases[net, f_?(needsInvSQ) :> momentumOf[f], Infinity];
  {env, nenv} = buildEnv[allMom, invMom, invSMom];

  ntLog["[prof] NumTrace analyseDiagram (", Length[diagrams], " diagrams): ",
    First@AbsoluteTiming[diags = analyseDiagram /@ diagrams], " s"];

  NTKernel[<|
    "Diagrams"  -> diags,
    "Env"       -> env,
    "NEnv"      -> nenv,
    "Frame"     -> frame,
    "Args"      -> args,
    "Dressings" -> dress
  |>]
];

(* One diagram -> {pure-scalar coeff, axis-id map, tensor components}. Components keep
   their factors un-expanded (heads, Plus-vertices, Times-structures); the recursive
   et compiler in Codegen turns Plus -> et::add, Times -> contract_all. *)
analyseDiagram[diagram_] := Module[{factors, tensorF, ids},
  factors = rewriteDressedNums @ splitSelfTraces[If[Head[diagram] === Times, List @@ diagram, {diagram}]];
  tensorF = Select[factors, ! scalarQ[#] &];
  (* Partition labels by sector: spinor (Dirac) axes get a disjoint high id range (>=100)
     so the et engine — which contracts axes by matching id — never fuses a spinor axis
     with a Lorentz/colour axis that happened to be numbered the same. *)
  ids     = With[{labs = DeleteDuplicates @ Flatten[allLabels /@ tensorF],
                  spn  = DeleteDuplicates @ Flatten[allSpinorLabels /@ tensorF]},
              With[{nonsp = DeleteCases[labs, Alternatives @@ spn]},
                Join[AssociationThread[nonsp -> Range[0, Length[nonsp] - 1]],
                     AssociationThread[spn -> Range[100, 99 + Length[spn]]]]]];
  <|
    (* contractFlavour collapses any flavour-δ chain that the dressing collection factored out of an
       eager numerator into the coeff (it straddled the Plus before, so FromFunKit could not). A no-op
       when collection is off / no flavour δ survives. *)
    "Coeff"      -> If[TrueQ[$ntDressCollect], contractFlavour[Times @@ Select[factors, scalarQ]],
                       Times @@ Select[factors, scalarQ]],
    "Ids"        -> ids,
    (* "Constant" is consumed by Codegen's per-component dispatch, where it means not merely
       "momentum-free" but "a constant SU(N) component": it is handed to compileColG, which only
       understands group heads. A momentum-free DIRAC structure (a bare closed spinor δ-loop, tr[1] = 4,
       with no ntVec/projector to carry momentum) is momentum-free yet carries no colour meaning, so
       compileColG emitted it as raw Mathematica (`colFacG[ntDeltaDirac[d2,d3], <|a1 -> 0, ...|>]`)
       into the generator .cpp — an instant clang failure. Reproduce with <P_2,T_2> of AqbqDirect8.
       Listing the Dirac heads here makes such a component NON-constant, routing it to
       splitColourGroupsInv/compileDirac like any other Dirac structure. NOTE this is only correct
       together with the collapsed-loop tr(1)=4 restoration in compileDirac (see nEmptyLoops there):
       orderDiracFacs drops δ connectors, so the loop arrives token-free and the runtime's split_loops
       would otherwise discard it, giving a kernel 4x too small. *)
    "Components" -> (<|"Factors" -> orderFactors[#],
                       "Constant" -> FreeQ[#, _ntVec | _ntTransProj | _ntLongProj |
                                              _ntElectricProj | _ntMagneticProj | _ntDressedNum |
                                              _ntGamma | _ntGamma5 | _ntSigma | _ntDeltaDirac]|> &
                     /@ connectedComponents[tensorF])
  |>
];

(* ---- per-diagram label validation ------------------------------------------------
   Run on ONE diagram AFTER expandBridges (a flat Times whose factors are heads, eager
   single-sector Plus vertices, Powers and scalars). The invariant, per diagram:
     * a label occurring ONCE  is a free (external) index of the whole trace;
     * a label occurring TWICE is contracted — the et engine pairs the two axes by id;
     * a label occurring 3+ times is MALFORMED. analyseDiagram assigns ONE id per DISTINCT
       label, make_eplan (axplan.hpp) pairs axes first-match-wins between exactly TWO
       operands, and contract_all is a LEFT FOLD — so four axes sharing an id silently
       contract as two independent pairs. No crash, just a wrong number. That is why this
       must Abort rather than Message.
   A Plus factor needs care and is why the old whole-net Cases[...,Infinity] count was
   useless: an eager vertex sum's summands LEGITIMATELY repeat the same free labels. A Plus
   is counted ONCE, via the free set its summands must agree on — that agreement is itself
   the et::add alignment precondition freeIdx[_Plus] assumes without checking. A summand's
   internal dummies are private and must not appear anywhere else. *)

NumTrace::badlabel = "Diagram `1`: index label `2` occurs `3` times (expected 1 = free, \
2 = contracted). The et engine contracts axes by matching id, so `3` axes sharing this label \
are silently mis-paired into a wrong number. Offending diagram:\n`4`";
NumTrace::plusfree = "Diagram `1`: the summands of an eager (un-distributed) sum expose \
DIFFERENT free indices `2` — et::add cannot align them. Offending sum:\n`3`";
NumTrace::privclash = "Diagram `1`: label(s) `2` are private dummies of one factor but also \
occur outside it — a dummy-name collision between two independently generated objects. \
Offending diagram:\n`3`";

(* {exposed-multiset, private-set, bad-list} of a (sub)expression. *)
labelCensus[e_] := Which[
  tensorQ[e],
    With[{ls = labelsOf[e]},
      (* a label repeated WITHIN one head is a legal self-trace; splitSelfTraces resolves it *)
      {DeleteDuplicates[ls], Cases[Tally[ls], {l_, c_} /; c >= 2 :> l], {}}],

  Head[e] === Power && IntegerQ[e[[2]]] && e[[2]] >= 1 && ! scalarQ[e],
    (* n copies sharing the SAME labels: a closed self-contraction (see compileTInv) *)
    Module[{c = labelCensus[e[[1]]], n = e[[2]]},
      Which[
        n === 1, c,
        n === 2, {{}, Union[c[[1]], c[[2]]], c[[3]]},
        True,    {{}, Union[c[[1]], c[[2]]], Join[c[[3]], ({#, n} &) /@ c[[1]]]}]],

  Head[e] === Times,
    Module[{sub = labelCensus /@ Select[List @@ e, ! scalarQ[#] &], exp, tal, priv, bad},
      If[sub === {}, Return[{{}, {}, {}}]];
      exp  = Join @@ sub[[All, 1]];
      tal  = Tally[exp];
      priv = Union @@ sub[[All, 2]];
      bad  = Join[Join @@ sub[[All, 3]], Cases[tal, {l_, c_} /; c > 2 :> {l, c}]];
      (* a child's private dummy must not be exposed or private anywhere else *)
      Do[With[{mine = sub[[i, 2]],
               others = Union[Join @@ Delete[sub, i][[All, 1]], Join @@ Delete[sub, i][[All, 2]]]},
           With[{clash = Intersection[mine, others]},
             If[clash =!= {}, bad = Join[bad, ({#, "private-clash"} &) /@ clash]]]],
         {i, Length[sub]}];
      {Cases[tal, {l_, c_} /; c == 1 :> l],
       Union[priv, Cases[tal, {l_, c_} /; c == 2 :> l]],
       bad}],

  Head[e] === Plus,
    Module[{sub = labelCensus /@ (List @@ e), frees},
      frees = Sort /@ sub[[All, 1]];
      {If[frees === {}, {}, First[frees]],
       Union @@ sub[[All, 2]],
       Join[Join @@ sub[[All, 3]],
            If[Length[DeleteDuplicates[frees]] > 1, {{frees, "plus-free-mismatch"}}, {}]]}],

  True, {{}, {}, {}}];

(* Escape hatch: NT_NO_LABEL_CHECK=1 disables (the census is O(net), not a hot path). *)
$ntCheckLabels = (Environment["NT_NO_LABEL_CHECK"] === $Failed);

checkLabels[diagram_, idx_] := If[! TrueQ[$ntCheckLabels], {},
  Module[{c = labelCensus[diagram], bad},
    bad = c[[3]];
    If[bad =!= {},
      Do[Switch[b[[2]],
           "plus-free-mismatch", Message[NumTrace::plusfree, idx, b[[1]], Short[diagram, 8]],
           "private-clash",      Message[NumTrace::privclash, idx, b[[1]], Short[diagram, 8]],
           _,                    Message[NumTrace::badlabel, idx, b[[1]], b[[2]], Short[diagram, 8]]],
         {b, bad}];
      Abort[]];
    c[[1]]]];
