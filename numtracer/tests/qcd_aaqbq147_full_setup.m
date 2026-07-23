(* ::Package:: *)

(* qcd_aaqbq147_full_setup.m — TEST fixture. Sibling of qcd_aaqbq147_setup.m, differing in ONE thing:
   that file restricts the external basis to element 2 (TBRestrictBasis[..., "transAAqbqR", {2}]),
   while this one keeps the FULL production transAAqbqMinimal (elements 1 AND 2), so ZAAqbq1 projects
   element 1 and ZAAqbq2 element 2. This is the fixture the lever-(b) validation runs against — see
   NUMTRACER_LEVER_B_BRIEF.md — and the flow that was intractable before it (>8.66 CPU-hr, OOM at
   20 GB, still inside phase A).

   Both are the FULL-BASIS sibling of qcd_aaqbq_setup.m: identical except
   the INTERNAL quark-gluon vertex {A,qb,q} uses the full "AqbqDirect147" tensor basis {1,4,7} (with
   dressings ZAqbq1/4/7) instead of the struct-1-only "AqbqDirect1". This is the Stage-4 north-star
   configuration: with COLLECTION ON the three internal-vertex structures fold to ONE collected
   ntDiracSlot per leg rather than distributing into 3^n diagrams — the payoff the collected-slot work
   exists to enable. Everything else (fields, propagators, the transAAqbqR external vertex restricted
   to element 2 to keep the proxy small+dense, dressing rules, frames) matches qcd_aaqbq_setup.m.
   Cloned rather than folded in for the same reason qcd_aaqbq_setup.m / qcd147_setup.m were: a
   truncation/basis change perturbs every frozen kernel that shares the setup.

   WHY THIS FIXTURE EXISTS. It is the only thing in the test suite with real sub-term density. The
   generator contracts one Dirac trace per (net, sub-term); the whole point of the global sub-term
   dedup (Codegen.m) is that the same trace recurs across nets. Every OTHER flow in tests/gen has
   exactly ONE sub-term per net, so the dedup is a no-op there and the merge, the trace cache and the
   tree fold are all exercised by nothing. This one has hundreds per net, and real redundancy.

   WHERE THE DENSITY COMES FROM (it is not obvious, and three plausible guesses are wrong):
   sub-terms are the monomials of Expand[Π (Dirac-bearing Plus factors)] (Codegen.m,
   splitColourGroupsInv). A Dirac sum only reaches that Expand if it survives three filters:
     - it is NOT a dressed structure sum (those become separate DIAGRAMS, not sub-terms — so the
       SIZE OF A VERTEX'S TENSOR BASIS DOES NOT DRIVE DENSITY: AqbqDirect147 has three structures
       and still yields 1 sub-term/net);
     - it is NOT a collectible propagator numerator (those collapse into one dtslot);
     - it is NOT a foldable 2-term gamma commutator (foldDiracSigma collapses those to one sigma
       token — which is exactly why the existing flows stay at 1).
   What survives all three is a multi-gamma structure with NUMERIC internal coefficients that is not
   a commutator. In this basis set that is precisely Tmunurho = {[g_mu,g_nu], g_rho}, the building
   block of transAAqbqMinimal's element 2 (T^(R)) — and it enters the trace TWICE, once as the
   in-loop AAqbq vertex and once inside the projector's dual. That product is the density.

   Hence the restriction below keeps element {2} and drops element {1}: restricting to {1} instead
   would DESTROY the density (T^(B) is a bare commutator Tmunu = [g_mu,g_nu], which sigma-folds to a
   single token). Keeping only {2} also makes the Gram 1x1, which is what keeps this cheap. *)

(* Pin FunKit's backend BEFORE any derivative/trace is taken: FRoute is leg-order sensitive, so
   the backend choice silently changes the loop-momentum routing and thus every frozen test
   integrand. See backend_pin.m. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "backend_pin.m"}]];

fields = <|
  "Commuting" -> {A[p, {v, c}]},
  "Grassmann" -> {{cb[p, {c}], c[p, {c}]}, {qb[p, {d, A, F}], q[p, {d, A, F}]}}
|>;

(* FULL internal quark-gluon basis {1,4,7} (this is the AqbqDirect147 switch); AqbqDirect1 kept
   registered too, matching qcd147_setup.m. *)
TBRestrictBasis["AqbqDirect", "AqbqDirect147", {1, 4, 7}];
TBRestrictBasis["AqbqDirect", "AqbqDirect1", {1}];
(* FULL transAAqbqMinimal (elements 1 AND 2) — the production ZAAqbq basis, so ZAAqbq1 projects
   element 1 (the light T^(B) bare commutator). No element restriction. *)

truncation = <|
  GammaN -> {{A, A}, {A, A, A}, {A, A, A, A}, {A, cb, c}, {cb, c}, {A, qb, q}, {A, A, qb, q}, {qb, q}},
  Propagator -> {{A, A}, {cb, c}, {qb, q}},
  R -> {{A, A}, {cb, c}, {qb, q}},
  Rdot -> {{A, A}, {cb, c}, {qb, q}},
  Field -> {{}}
|>;

bases = <|
  GammaN -> {{A, A} -> {"AA", 1}, {A, A, A} -> "AAAClass", {A, A, A, A} -> "AAAAClass",
     {A, cb, c} -> {"Acbc", 1}, {cb, c} -> "cbc", {A, qb, q} -> "AqbqDirect147",
     {A, A, qb, q} -> "transAAqbqMinimal", {qb, q} :> "qbq"},
  Propagator -> {{A, A} -> {"AA", 1}, {cb, c} -> "cbc", {qb, q} -> "qbq"},
  R -> {{A, A} -> {"AA", 1}, {cb, c} -> "cbc", {qb, q} -> {"qbq", 1}},
  Rdot -> {{A, A} -> {"AA", 1}, {cb, c} -> "cbc", {qb, q} -> {"qbq", 1}}
|>;

Setup = <|
  "FieldSpace" -> fields,
  "Truncation" -> truncation,
  "FeynmanRules" -> bases,
  "DiagramStyling" -> <|"Styles" -> {A -> {Orange}, c -> {Black, Dashed}, q -> {Black}}|>
|>;
FSetGlobalSetup[Setup];

SP3Patt[p1e_, p2e_, p3e_] := {Sqrt[(sp[p1, p1] + sp[p2, p2] + sp[p3, p3])/3]} /.
    {p1 :> p1e, p2 :> p2e, p3 :> p3e} // UseLorentzLinearity // FullSimplify;
SP4Patt[p1e_, p2e_, p3e_, p4e_] := {Sqrt[(sp[p1, p1] + sp[p2, p2] + sp[p3, p3] + sp[p4, p4])/4]} /.
    {p1 :> p1e, p2 :> p2e, p3 :> p3e, p4 :> p4e} // UseLorentzLinearity // FullSimplify;

dressingRules[expr_] := Module[{tmp},
  tmp = ReplaceRepeated[expr, {
    (* Propagators *)
    dressing[GammaN, {cb, c}, 1, {p1_, p2_}] :> -Zc[Sqrt[sp[p2, p2]]] sp[p2, p2],
    dressing[GammaN, {A, A}, 1, {p1_, p2_}] :> ZA[Sqrt[sp[p2, p2]]] sp[p2, p2],
    dressing[GammaN, {qb, q}, 1, {p1_, p2_}] :> -Zq[Sqrt[sp[p2, p2]]],
    dressing[GammaN, {qb, q}, 2, {p1_, p2_}] :> -Mq[Sqrt[sp[p2, p2]]],
    (* Inverse propagators (contain the regulator) *)
    dressing[InverseProp, {cb, c}, 1, {p1_, p2_}] :> -(Zc[Sqrt[sp[p2, p2]]] sp[p2, p2] + RB[k^2, sp[p2, p2]] Zc[k]),
    dressing[InverseProp, {A, A}, 1, {p1_, p2_}] :> ZA[Sqrt[sp[p2, p2]]] sp[p2, p2] + RB[k^2, sp[p2, p2]] ZA[evP],
    dressing[InverseProp, {qb, q}, 1, {p1_, p2_}] :> -(Zq[Sqrt[sp[p2, p2]]] + (Zq[k] RF[k^2, sp[p2, p2]])/Sqrt[sp[p2, p2]]),
    dressing[InverseProp, {qb, q}, 2, {p1_, p2_}] :> -Mq[Sqrt[sp[p2, p2]]],
    (* Strong couplings *)
    dressing[GammaN, {A, cb, c}, 1, {p1_, p2_, p3_}] :> ZAcbc[p1, p2],
    dressing[GammaN, {A, A, A}, 1, {p1_, p2_, p3_}] :> ZA3[p1, p2],
    dressing[GammaN, {A, A, A, A}, 1, {p1_, p2_, p3_, p4_}] :> ZA4[p1, p2, p3],
    dressing[GammaN, {A, qb, q}, 1, {p1_, p2_, p3_}] :> ZAqbq1[p1, p2],
    dressing[GammaN, {A, qb, q}, 2, {p1_, p2_, p3_}] :> ZAqbq4[p1, p2],
    dressing[GammaN, {A, qb, q}, 3, {p1_, p2_, p3_}] :> ZAqbq7[p1, p2],
    (* the dense vertex: basis restricted to element 2 (T^(R)), so it carries dressing index 1 *)
    dressing[GammaN, {A, A, qb, q}, 1, {p1_, p2_, p3_, p4_}] :> ZAAqbq1[p1, p2, p3],
    dressing[GammaN, {A, A, qb, q}, 2, {p1_, p2_, p3_, p4_}] :> ZAAqbq2[p1, p2, p3],
    (* Parametrizations *)
    ZAcbc[p1_, p2_] :> ZAcbc @@ SP3Patt[p1, p2, -p1 - p2],
    ZA3[p1_, p2_] :> ZA3 @@ SP3Patt[p1, p2, -p1 - p2],
    ZA4[p1_, p2_, p3_] :> ZA4 @@ SP4Patt[p1, p2, p3, -p1 - p2 - p3],
    ZAqbq1[p1_, p2_] :> ZAqbq1 @@ SP3Patt[p1, p2, -p1 - p2],
    ZAqbq4[p1_, p2_] :> ZAqbq4 @@ SP3Patt[p1, p2, -p1 - p2],
    ZAqbq7[p1_, p2_] :> ZAqbq7 @@ SP3Patt[p1, p2, -p1 - p2],
    ZAAqbq1[p1_, p2_, p3_] :> ZAAqbq1 @@ SP4Patt[p1, p2, p3, -p1 - p2 - p3],
    ZAAqbq2[p1_, p2_, p3_] :> ZAAqbq2 @@ SP4Patt[p1, p2, p3, -p1 - p2 - p3],
    nZA -> 6,
    evP :> (k^nZA + 1)^(1/nZA),
    devP :> k^(-1 + nZA) (1 + k^nZA)^(-1 + 1/nZA),
    (* Regulator derivatives *)
    dressing[Rdot, {A, A}, 1, {p1_, p2_}] :> ZA[evP] RBdot[k^2, sp[p2, p2]] +
       RB[k^2, sp[p2, p2]] (dtZA[evP] + k*devP*(ZA[1.02 evP] - ZA[evP])/(0.02*evP)),
    dressing[Rdot, {cb, c}, 1, {p1_, p2_}] :> Zc[k] RBdot[k^2, sp[p2, p2]] +
       RB[k^2, sp[p2, p2]] (dtZc[k] + k (Zc[1.02*k] - Zc[k])/(0.02*k)),
    (* Minus sign here due to the vertex basis: it's pdash[p1,d1,d2]! *)
    dressing[Rdot, {qb, q}, 1, {p1_, p2_}] :> -(Zq[k] RFdot[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]] +
       RF[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]] (dtZq[k] + k (Zq[1.02*k] - Zq[k])/(0.02*k))),
    (* Regulators *)
    dressing[R, {qb, q}, 1, {p1_, p2_}] :> -(Zq[k] RF[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]])
  }];
  UseLorentzLinearity[tmp] //. {lf1 -> l1} // Return;
];

FSetSymmetricDressing[GammaN, {A, A}];
FSetSymmetricDressing[GammaN, {cb, c}];
FSetSymmetricDressing[GammaN, {qb, q}];

(* vacuum multi-angle scalar reduction for 3-/4-point vertices: sp/cos -> {l1, p, cosl1p_i}. *)
SPParam[expr_] := UseLorentzLinearity[expr] //. {
  lf1 -> l1,
  sp[p, p] -> p^2, sp[l1, l1] -> l1^2,
  sp[l1, p1] -> p l1 cos[l1, p1], sp[l1, p2] -> p l1 cos[l1, p2],
  sp[l1, p3] -> p l1 cos[l1, p3], sp[l1, p4] -> p l1 cos[l1, p4],
  Sqrt[a_^2] :> a, (a_^2)^(n_/2) :> a^n, Power[Power[l1_, 2], Rational[n_, 2]] :> l1^n,
  cos[l1, p1] :> cosl1p1, cos[l1, p2] :> cosl1p2,
  cos[l1, p3] :> cosl1p3, cos[l1, p4] :> cosl1p4
};

(* Compile-time SU(N) sizes: colour SU(3), flavour SU(2). transAAqbqMinimal requires BOTH groups. *)
SetNc[3];
SetNf[2];
$Assumptions = k > 0 && p > 0 && l1 > 0 && -1 < cos1 < 1 && -1 < cos2 < 1;

$ntShimOpts = {"RuntimeInclude" -> "shim.hpp", "ExtraIncludes" -> {},
  "KernelNamespace" -> "DiFfRG", "SupportNamespace" -> "DiFfRG",
  "DressingType" -> "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>"};
