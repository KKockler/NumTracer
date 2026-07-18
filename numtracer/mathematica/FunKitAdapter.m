(* ::Package:: *)

(* FunKit adapter: rewrite a traced flow expression (FunKit's FormTracer head
   vocabulary) into the NumTracer DSL. A pure, local rewrite over a small closed
   set of heads — it does NOT touch the TensorBases internals.

   Expected input is a flow AFTER `// dressingRules // PropParam` (or SPParam): the
   dressings are ZA/Zc/RB expressions and the scalar dot products are already
   reduced to the runtime scalars (l1, p, cos1, ...). So only the TENSOR heads are
   rewritten here; their momenta are resolved by the frame in MakeKernel. (`sps`,
   the finite-T scalar product, is its own symbol and left untouched — finite-T
   support comes later.)

       flow = traceExprcbc // dressingRules // PropParam;
       net  = FromFunKit[flow];

   Structure:  FEx[...] = sum of terms,  FTerm[...] = product of factors.

   We dispatch on the head's *name* (SymbolName), not the symbol itself: the FunKit
   heads live in TensorBases`/FunKit` contexts not necessarily on $ContextPath when
   this file loads, so a literal `transProj` pattern would be a different symbol. *)

(* The FunKit-token -> NumTracer-head map, EXCEPT the SU(N) group tokens: those carry no rank
   of their own, so FromFunKit injects the rank N (a bound integer) into the new N-parameterized
   heads — see $sunMap below. This part is rank-independent. *)
$ffMap = <|
  "FEx" -> Plus, "FTerm" -> Times, "NonCommutativeMultiply" -> Times,
  "transProj" -> ntTransProj, "longProj" -> ntLongProj, "vec" -> ntVec,
  "deltaLorentz" -> ntMetric,
  (* Dirac (spinor) sector. FunKit emits a slashed momentum as gamma[mu,..] * vec[q,mu]
     (the gamma carries the Lorentz axis), so ntGamma -> Dirac::gamma_axis. *)
  "gamma" -> ntGamma, "gamma5" -> ntGamma5, "deltaDirac" -> ntDeltaDirac,
  (* flavour-TRIVIAL Kronecker delta -> a private head, contracted to a power of Nf below (its
     dimension Nf is symbolic). Correct ONLY for a flavour-blind closed loop (delta^{ii} = Nf,
     e.g. Zq/ZA quark loops). A fundamental flavour delta sitting INSIDE a τ-trace must instead
     be emitted as deltaFlavFundGen (-> ntSUNT/ntSUNDeltaFund via $sunMap) so it stays in the
     engine and does not break the trace by collapsing. *)
  "deltaFundFlav" -> flavDelta,
  (* any scalar product the notebook's reduction leaves (e.g. external-external SP
     constants sp[p_i,p_j] that SPParam does not touch) -> ntSP, resolved by the frame
     as the Euclidean dot of the components. For Zc/ZA nothing survives PropParam, so
     this is a no-op there. *)
  "sp" -> ntSP,
  (* finite-T spatial scalar product. Maps to ntSPS, resolved by the frame as the SPATIAL dot
     (components 1..3) of the momenta. FunKit's component access vec[q, 0] (a literal integer index)
     rides the "vec" -> ntVec map above and is routed to the scalar temporal component q_0 by the
     integer-index classification in DSL.m — no special rule needed here. *)
  "sps" -> ntSPS
|>;

(* The SU(N) group tokens, mapped to the N-parameterized heads with the group rank `n` injected
   as the leading argument. `nc` is the colour rank Global`Nc; `nf` the isospin rank (FromFunKit's
   "FlavourGroup" option). Colour structure (FCol/deltaAdjCol/TCol/deltaFundCol) builds against
   nc; the hand-rolled QM-model isospin tokens (the τ Yukawa generator, the pion f^{abc}
   self-coupling, the adjoint/fundamental isospin deltas) build against nf. Both go through the
   SAME four heads — the engine separates the groups by their disjoint contraction ids. *)
sunMap[nc_, nf_] := <|
  "FCol" -> (ntSUNf[nc, ##] &), "deltaAdjCol" -> (ntSUNDeltaAdj[nc, ##] &),
  "TCol" -> (ntSUNT[nc, ##] &), "deltaFundCol" -> (ntSUNDeltaFund[nc, ##] &),
  "fFlav" -> (ntSUNf[nf, ##] &), "deltaAdjFlav" -> (ntSUNDeltaAdj[nf, ##] &),
  "tauFlav" -> (ntSUNT[nf, ##] &), "deltaFlavFundGen" -> (ntSUNDeltaFund[nf, ##] &),
  (* FUNDAMENTAL Levi-Civita: N indices (colour SU(3) -> 3, isospin SU(2) -> 2). Contracted into
     Kronecker deltas by expandFundEps in DSL.m — no engine token. *)
  "epsFundCol" -> (ntEpsFund[nc, ##] &), "epsFundFlav" -> (ntEpsFund[nf, ##] &),
  (* ADJOINT Levi-Civita, SU(2) ONLY — see adjEps. *)
  "epsAdjCol" -> (adjEps[nc, ##] &), "epsAdjFlav" -> (adjEps[nf, ##] &)
|>;

(* ---- adjoint Levi-Civita: SU(2) ONLY ---------------------------------------------------------
   The adjoint of SU(N) has dimension N^2-1, so its epsilon carries N^2-1 indices: 3 for SU(2), but
   8 for SU(3). ShowFormTracerDefinitions[] displays `epsAdjCol[a, b, c]` with a generic THREE-index
   signature, which makes a 3-index adjoint epsilon at Nc=3 look legal. It is not — that is the
   likely mistake, so the refusal below names it.

   !! SU(2) ONLY !!  At rank 2 (and ONLY at rank 2) the adjoint epsilon coincides with the structure
   constant: eps^{abc} = f^{abc}, exactly, coefficient +1. This is NOT a general fact — for SU(3),
   f^{abc} is not an epsilon at all (it is not totally antisymmetric in the same sense and has
   different nonzero entries), so nothing here may be read as an SU(N) statement.

   The +1 is derived from NumTracer's OWN conventions, not assumed: sun_net.hpp:233-247 builds the
   N=2 generators as T^a = sigma^a/2 in (x,y,z) order, and sun_net.hpp:225 defines
   f^{abc} = -2i tr([T^a,T^b] T^c). With [T^a,T^b] = i eps^{abc} T^c and tr(T^c T^d) = delta^{cd}/2
   this gives f^{abd} = eps^{abd}. It is still pinned by a test (a contraction LINEAR in the
   coefficient, so a sign flip cannot hide) rather than trusted.

   Rewriting to ntSUNf reuses an already-tested primitive and — unlike a pair-contraction scheme —
   accepts an UNPAIRED adjoint epsilon (eps^{abc} T^b T^c and friends), which is common and
   perfectly well-defined. *)
(* FunKit's DECLARED head vocabulary (ShowFormTracerDefinitions[]). The closed-world list the guard
   in FromFunKit checks against: every one of these must be mapped, or refused on purpose. *)
$funKitHeads = {"FEx", "FTerm", "deltaLorentz", "vec", "sp", "sps",
  "deltaDirac", "gamma", "gamma5", "sigma", "transProj", "longProj",
  "deltaAdjCol", "deltaFundCol", "FCol", "TCol", "epsAdjCol", "epsFundCol",
  "deltaAdjFlav", "deltaFundFlav", "fFlav", "tauFlav", "TFlav",
  "epsAdjFlav", "epsFundFlav", "deltaFlavFundGen", "epsLorentz"};
(* Handled outside `map`: TFlav by the hasIso rewrite above. *)
$ffHandledElsewhere = {"TFlav"};
FromFunKit::untranslated = "the FunKit head(s) `1` appear in the input but have no entry in \
$ffMap / sunMap. An untranslated head does NOT fail loudly downstream: DSL.m's scalarQ is a FreeQ \
over the KNOWN nt* heads, so an unknown head is classified as a SCALAR COEFFICIENT — its indices \
become invisible to labelsOf/freeIdx, the diagram reports spurious free (open) legs, checkLabels \
accepts it (open legs are legal), and the raw head is CForm'd into the generated C++. This is how \
epsFundCol/epsFundFlav went undetected. Add a $ffMap/sunMap entry, or refuse the input explicitly.";
(* epsLorentz is REFUSED on purpose, not mapped: it is the 3D SPATIAL epsilon (O(3) after the
   heat-bath split), whereas ntEpsilon is 4D and hard-wired to four labels (DSL.m labelsOf,
   Codegen.m builderInv and the Length[lst] == 4 reconstruction). Mapping one onto the other is a
   silent dimension error. The pair route is cheap (2 terms at D=3) but produces SPATIAL deltas, and
   NumTracer has no spatial-delta head — ntMetric is the 4D Euclidean metric and would wrongly
   include the temporal component. So the blocker is a missing head, not cost: the fix is a spatial
   delta, not more epsilon machinery. It falls through to FromFunKit::untranslated. *)

FromFunKit::epsadj = "Adjoint Levi-Civita at SU(`1`) with `2` indices. NumTracer supports the \
adjoint epsilon ONLY at rank 2, where eps^abc coincides exactly with the structure constant f^abc \
(T^a = sigma^a/2, f = -2i tr([T^a,T^b]T^c); see sun_net.hpp:225) and is rewritten to ntSUNf[2,...]. \
That identification is SU(2)-SPECIFIC and does NOT generalise. At SU(`1`) the adjoint epsilon \
carries `3` indices, not `2` — ShowFormTracerDefinitions[] shows epsAdjCol[a,b,c] only as a generic \
three-index illustration, and taking that arity literally at rank > 2 is the likely mistake here. \
Aborting rather than guessing: a rank/arity mismatch would silently contract against the wrong \
group.";
adjEps[n_, idx__] := If[n === 2 && Length[{idx}] === 3, ntSUNf[2, idx],
  Message[FromFunKit::epsadj, n, Length[{idx}], n^2 - 1]; Abort[]];

(* Contract the flavour Kronecker deltas: their indices are disjoint from every tensor
   sector, so a chain collapses (delta[x,y] delta[y,z] -> delta[x,z]) and a closed loop
   delta[x,x] -> Nf. The result is a scalar power of Nf that the per-diagram coefficient
   carries (cancelling the projector's 1/Nf for a flavour-trivial flow like Zq). *)
contractFlavour[e_] := e //. {
  flavDelta[x_, y_] flavDelta[y_, z_] :> flavDelta[x, z],
  flavDelta[x_, x_] :> Global`Nf,
  Power[flavDelta[x_, y_], 2] :> Global`Nf
};

(* "FlavourGroup" -> Automatic resolves the isospin SU(N) rank to Global`Nf when that is a defined
   integer, else defaults to 2; an explicit integer overrides it. The colour rank is Global`Nc. *)
Options[FromFunKit] = {"FlavourGroup" -> Automatic, "DressingCollection" -> True};

(* Rewrite heads (injecting each SU(N) group's rank N into the new heads), then DISTRIBUTE
   (expandBridges turns the propagator-numerator structure sums into separate flat-product
   diagrams), THEN contract the flavour deltas — only once each diagram is a flat product do a
   flavour chain's links sit in one Times so they can collapse to a power of Nf (before
   distribution they straddle a Plus and cannot).
   "DressingCollection" -> True sets the gate BEFORE expandBridges (FromFunKit runs it before
   NumTrace), so dressed Dirac numerators are kept eager here too; pass the SAME value to NumTrace. *)
FromFunKit[expr_, OptionsPattern[]] := Module[{nf, map, hasIso, isoRewritten, res},
  nf  = OptionValue["FlavourGroup"] /. Automatic :> If[IntegerQ[Global`Nf], Global`Nf, 2];
  map = Join[$ffMap, sunMap[Global`Nc, nf]];
  (* ISOSPIN GENERATORS (quark-meson flows). The notebook auxiliary `TFlav` is the SU(nf) FUNDAMENTAL
     flavour generator: TFlav[a, f1, f2] = (T^a)_{f1 f2} for an adjoint index a != 0, and the flavour
     SINGLET TFlav[0,f1,f2] = deltaFundFlav[f1,f2]/Sqrt[2 Nf]. It is NOT in $ffMap (it is a notebook
     symbol, not a FunKit token), so without this it leaks through as an opaque scalar and the isospin
     trace never closes — leaving the external flavour indices dangling (the pion (-I) factors then fail
     to cancel and a spurious imaginary part survives). Route the WHOLE fundamental-flavour sector into
     the SU(nf) engine: map TFlav to ntSUNT / ntSUNDeltaFund, and switch the connecting fundamental
     deltas from the blind `flavDelta` (which only collapses a genuinely CLOSED line to Nf) to the
     in-engine ntSUNDeltaFund so the generator trace tr(T^a ... T^a) actually contracts. Gated on the
     presence of TFlav, so flavour-blind flows (Zq/ZA/ZAqbq1/4/7, ...) are byte-identical. *)
  hasIso = ! FreeQ[expr, Global`TFlav];
  If[hasIso, map["deltaFundFlav"] = (ntSUNDeltaFund[nf, ##] &)];
  isoRewritten = If[hasIso,
    expr //. {Global`TFlav[0, f1_, f2_] :> ntSUNDeltaFund[nf, f1, f2]/Sqrt[2 Global`Nf],
              Global`TFlav[a_, f1_, f2_]  :> ntSUNT[nf, a, f1, f2]},
    expr];
  (* ---- THE ROOT-CLASS GUARD ------------------------------------------------------------------
     An untranslated FunKit head does NOT fail loudly downstream. DSL.m's scalarQ is a FreeQ over the
     KNOWN nt* heads, so an unknown head is classified as a SCALAR COEFFICIENT: its indices become
     invisible to labelsOf/freeIdx, the diagram reports spurious free (open) legs, checkLabels
     accepts it (open legs are legal), and the raw Mathematica head is CForm'd into the generated
     C++. That is exactly how epsFundCol/epsFundFlav went undetected through a whole debugging
     session — the visible symptom was 8 dangling indices three layers away from the cause.
     So: refuse any head from FunKit's declared vocabulary that has no entry in the map.
     NOTE the map is completed CONDITIONALLY just above (the hasIso branch promotes deltaFundFlav),
     so this must read `map`, never $ffMap — reading the wrong one is itself a way to be misled. *)
  With[{present = DeleteDuplicates @ Cases[isoRewritten, (h_Symbol)[___] :> SymbolName[h], {0, Infinity}]},
    With[{leftover = Complement[Intersection[present, $funKitHeads], Keys[map], $ffHandledElsewhere]},
      If[leftover =!= {}, Message[FromFunKit::untranslated, leftover]; Abort[]]]];
  $ntDressCollect = TrueQ[OptionValue["DressingCollection"]];
  ntLog["[prof] FromFunKit (head rewrite + expandBridges): ",
    First@AbsoluteTiming[res = contractFlavour @ expandBridges[
      isoRewritten //. (h_Symbol)[a___] /; KeyExistsQ[map, SymbolName[h]] :> map[SymbolName[h]][a]]], " s"];
  res];
