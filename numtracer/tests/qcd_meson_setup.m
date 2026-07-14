(* ::Package:: *)

(* qcd_meson_setup.m — TEST fixture (sibling of qcd_setup.m): the QCD + light-meson field space
   (gluon + ghost + quark + the pion ΠL and sigma σL of a quark-meson sector), its truncation,
   tensor bases, meson vertices, diagrammatic rules and dressing rules. Ported faithfully from
   ~/Documents/Uni/Code/qcd-codes/vacuum/QCD_NumTracer/basic/QCD.nb (cells 2-7), the notebook that
   derives the etaPiL / etaSigL pion/sigma anomalous dimensions. This is NOT part of the NumTracer
   frontend — it is the per-theory setup a generator script Gets so it can build the (un-traced)
   pion two-point loop. The Dirac/colour/flavour trace itself is taken downstream, by either the
   NumTracer numeric backend or DiFfRG's FormTracer.

   Compile-time SU(N): colour SU(3), flavour SU(2) (Nf=2). A closed quark loop folds to the integer
   Nf=2, so the generated kernels carry NO Nf parameter — both the numeric and the FORM kernel bake
   Nf=2 here and thus share one signature, which is what makes the cross-check a pointwise compare. *)

(* ---- field space + truncation (QCD.nb cell 2) -------------------------------------------------- *)
(* Pin FunKit's backend BEFORE any derivative/trace is taken: FRoute is leg-order sensitive, so
   the backend choice silently changes the loop-momentum routing and thus every frozen test
   integrand. See backend_pin.m. *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "backend_pin.m"}]];

fields = <|
  "Commuting" -> {A[p, {v, c}], ΠL[p, {f}], σL[p]},
  "Grassmann" -> {{cb[p, {c}], c[p, {c}]}, {qb[p, {d, A, F}], q[p, {d, A, F}]}}
|>;

TBRestrictBasis["AqbqDirect", "AqbqDirect147", {1, 4, 7}];
TBRestrictBasis["AqbqDirect", "AqbqDirect1", {1}];

truncation = <|
  GammaN -> {{A, A}, {A, A, A}, {A, A, A, A}, {A, cb, c}, {cb, c}, {A, qb, q}, {qb, q},
     {qb, q, ΠL}, {qb, q, σL}, {ΠL, ΠL}, {σL, σL}},
  Propagator -> {{A, A}, {cb, c}, {qb, q}, {ΠL, ΠL}, {σL, σL}},
  R -> {{A, A}, {cb, c}, {qb, q}, {ΠL, ΠL}, {σL, σL}},
  Rdot -> {{A, A}, {cb, c}, {qb, q}, {ΠL, ΠL}, {σL, σL}},
  Field -> {{σL}}
|>;

bases = <|
  GammaN -> {{A, A} -> {"AA", 1}, {A, A, A} -> "AAAClass", {A, A, A, A} -> "AAAAClass",
     {A, cb, c} -> {"Acbc", 1}, {cb, c} -> "cbc", {A, qb, q} -> "AqbqDirect1", {qb, q} :> "qbq"},
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

(* ---- meson vertices (QCD.nb cell 3) ------------------------------------------------------------
   The (pointlike) light-meson self-interactions, obtained by Taylor-expanding the effective
   potential VL about the symmetric point and functionally differentiating in the pion direction.
   For the pion two-point function only λmΠnσ (the pure-meson σ/π vertices) and the Yukawa rules in
   diagRules are needed; the resulting meson loops are momentum-independent and drop out of the p^2
   (wave-function) extraction, leaving the quark loop. *)
TFlav[0, f1_, f2_] := deltaFundFlav[f1, f2]/Sqrt[2*Nf];

FunctionalD[expr_, v : ((f_)[_] | {(f_)[_], _Integer}).., OptionsPattern[]] :=
  Internal`InheritedBlock[{f},
    f /: D[f[x_], f[y_], NonConstants -> {f}] := DiracDelta[x - y];
    f /: D[f, f[y_], NonConstants -> {f}] := DiracDelta[#1 - y] &;
    D[expr, v, NonConstants -> {f}]];

Taylor[func_, var_, order_] := func + Total[Table[
   ((1/n!)*Symbol["d" ~~ ToString[n] ~~ ToString[func]])*Product[var[], {i, 1, n}], {n, 1, order}]];

λmΠnσ[pions_, sigmas_] := Module[{action, sigmaDeriv},
  action = Taylor[VL, (1/2)*((sigma^2 - sigma0^2) + Module[{b}, Π[b]^2]) &, Length[pions] + sigmas];
  sigmaDeriv = D[action, {sigma, sigmas}];
  Simplify[If[Length[pions] >= 1,
     FunctionalD[sigmaDeriv, (Π[#1] &) /@ pions /. List -> Sequence] //.
       DiracDelta[(x_) - (y_)] -> deltaAdjFlav[x, y] //. DiracDelta[-1 + (y_)] -> deltaAdjFlav[1, y],
     sigmaDeriv] //. {Π[a_] -> 0, sigma0 -> sigma, sigma -> Sqrt[2*rhoL],
     Derivative[n_][V][a_] :> Symbol[StringJoin["d", ToString[n], "V"]], V[a_] -> V}]];

(* ---- diagrammatic rules (QCD.nb cell 4) -------------------------------------------------------- *)
diagRules = Join[FMakeDiagrammaticRules[], {
  Propagator[{ΠL, ΠL}, {{p1_, {f1_}}, {p2_, {f2_}}}] :>
    deltaAdjFlav[f1, f2]*(1/(zPiL[Sqrt[sp[p2, p2]]]*sp[p2, p2] + d1V + RB[k^2, sp[p2, p2]])),
  Propagator[{σL, σL}, {{p1_}, {p2_}}] :>
    1/(zSigL[Sqrt[sp[p2, p2]]]*sp[p2, p2] + d1V + 2*rhoL*d2V + RB[k^2, sp[p2, p2]]),
  Rdot[{ΠL, ΠL}, {{p1_, {f1_}}, {p2_, {f2_}}}] :>
    deltaAdjFlav[f1, f2]*(RBdot[k^2, sp[p2, p2]] - etaPiL*RB[k^2, sp[p2, p2]]),
  Rdot[{σL, σL}, {{p1_}, {p2_}}] :> RBdot[k^2, sp[p2, p2]] - etaSigL*RB[k^2, sp[p2, p2]],
  GammaN[{σL, σL, σL, σL}, {{p1_}, {p2_}, {p3_}, {p4_}}] :> λmΠnσ[{}, 4],
  GammaN[{σL, σL, σL}, {{p1_}, {p2_}, {p3_}}] :> λmΠnσ[{}, 3],
  GammaN[{σL, ΠL, ΠL}, {{p1_}, {p2_, {f2_}}, {p3_, {f3_}}}] :> λmΠnσ[{f2, f3}, 1],
  GammaN[{σL, σL, ΠL, ΠL}, {{p1_}, {p2_}, {p3_, {f3_}}, {p4_, {f4_}}}] :> λmΠnσ[{f3, f4}, 2],
  GammaN[{ΠL, ΠL, ΠL, ΠL}, {{p1_, {f1_}}, {p2_, {f2_}}, {p3_, {f3_}}, {p4_, {f4_}}}] :>
    λmΠnσ[{f1, f2, f3, f4}, 0],
  GammaN[{σL, qb, q}, {{p1_}, {p2_, {d2_, A2_, F2_}}, {p3_, {d3_, A3_, F3_}}}] :>
    hSigL[Sqrt[(sp[p1, p1] + sp[p2, p2] + sp[p3, p3])/3]]*(deltaFundFlav[F2, F3]/Sqrt[2*Nf])*
      deltaDirac[d2, d3]*deltaFundCol[A2, A3],
  GammaN[{ΠL, qb, q}, {{p1_, {f1_}}, {p2_, {d2_, A2_, F2_}}, {p3_, {d3_, A3_, F3_}}}] :>
    hPiL[Sqrt[(sp[p1, p1] + sp[p2, p2] + sp[p3, p3])/3]]*TFlav[f1, F2, F3]*
      ((-I)*gamma5[d2, d3])*deltaFundCol[A2, A3]}];

(* ---- symmetric-point parametrisations + dressing rules (QCD.nb cell 5) ------------------------- *)
SP3Patt[p1e_, p2e_, p3e_] := FullSimplify[UseLorentzLinearity[
   {Sqrt[(sp[p1, p1] + sp[p2, p2] + sp[p3, p3])/3]} /. {p1 :> p1e, p2 :> p2e, p3 :> p3e}]];
SP4Patt[p1e_, p2e_, p3e_, p4e_] := FullSimplify[UseLorentzLinearity[
   {Sqrt[(sp[p1, p1] + sp[p2, p2] + sp[p3, p3] + sp[p4, p4])/4]} /. {p1 :> p1e, p2 :> p2e, p3 :> p3e, p4 :> p4e}]];

dressingRules[expr_] := Module[{tmp},
  tmp = expr //. {
    (* Propagators *)
    dressing[GammaN, {cb, c}, 1, {p1_, p2_}] :> (-Zc[Sqrt[sp[p2, p2]]])*sp[p2, p2],
    dressing[GammaN, {A, A}, 1, {p1_, p2_}] :> ZA[Sqrt[sp[p2, p2]]]*sp[p2, p2],
    dressing[GammaN, {qb, q}, 1, {p1_, p2_}] :> -Zq[Sqrt[sp[p2, p2]]],
    dressing[GammaN, {qb, q}, 2, {p1_, p2_}] :> -Mq[Sqrt[sp[p2, p2]]],
    (* Inverse propagators (contain the regulator) *)
    dressing[InverseProp, {cb, c}, 1, {p1_, p2_}] :> -(Zc[Sqrt[sp[p2, p2]]]*sp[p2, p2] + RB[k^2, sp[p2, p2]]*Zc[k]),
    dressing[InverseProp, {A, A}, 1, {p1_, p2_}] :> ZA[Sqrt[sp[p2, p2]]]*sp[p2, p2] + RB[k^2, sp[p2, p2]]*ZA[evP],
    dressing[InverseProp, {qb, q}, 1, {p1_, p2_}] :> -(Zq[Sqrt[sp[p2, p2]]] + (Zq[k]*RF[k^2, sp[p2, p2]])/Sqrt[sp[p2, p2]]),
    dressing[InverseProp, {qb, q}, 2, {p1_, p2_}] :> -Mq[Sqrt[sp[p2, p2]]],
    (* Strong couplings *)
    dressing[GammaN, {A, cb, c}, 1, {p1_, p2_, p3_}] :> ZAcbc[p1, p2],
    dressing[GammaN, {A, A, A}, 1, {p1_, p2_, p3_}] :> ZA3[p1, p2],
    dressing[GammaN, {A, A, A, A}, 1, {p1_, p2_, p3_, p4_}] :> ZA4[p1, p2, p3],
    dressing[GammaN, {A, qb, q}, 1, {p1_, p2_, p3_}] :> ZAqbq1[p1, p2],
    dressing[GammaN, {A, qb, q}, 2, {p1_, p2_, p3_}] :> ZAqbq4[p1, p2],
    dressing[GammaN, {A, qb, q}, 3, {p1_, p2_, p3_}] :> ZAqbq7[p1, p2],
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
    devP :> k^(-1 + nZA)*(1 + k^nZA)^(-1 + 1/nZA),
    (* Regulator derivatives *)
    dressing[Rdot, {A, A}, 1, {p1_, p2_}] :> ZA[evP]*RBdot[k^2, sp[p2, p2]] +
       RB[k^2, sp[p2, p2]]*(dtZA[evP] + k*devP*((ZA[1.02*evP] - ZA[evP])/(0.02*evP))),
    dressing[Rdot, {cb, c}, 1, {p1_, p2_}] :> Zc[k]*RBdot[k^2, sp[p2, p2]] +
       RB[k^2, sp[p2, p2]]*(dtZc[k] + k*((Zc[1.02*k] - Zc[k])/(0.02*k))),
    (* Minus sign here due to the vertex basis: it's pdash[p1,d1,d2]! *)
    dressing[Rdot, {qb, q}, 1, {p1_, p2_}] :> -(Zq[k]*(RFdot[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]]) +
       (RF[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]])*(dtZq[k] + k*((Zq[1.02*k] - Zq[k])/(0.02*k)))),
    (* Regulators *)
    dressing[R, {qb, q}, 1, {p1_, p2_}] :> -(Zq[k]*(RF[k^2, sp[p2, p2]]/Sqrt[sp[p2, p2]])),
    Mq[x_] :> hSigL[x]*(Sqrt[2*rhoL]/Sqrt[2*Nf])};
  UseLorentzLinearity[tmp] //. {lf1 -> l1}];

FSetSymmetricDressing[GammaN, {A, A}];
FSetSymmetricDressing[GammaN, {cb, c}];
FSetSymmetricDressing[GammaN, {qb, q}];

(* ---- scalar reductions (QCD.nb cell 6) --------------------------------------------------------- *)
(* vacuum 1-angle reduction for the propagators: sp/cos -> {l1, p, cos1} *)
PropParam[expr_] := UseLorentzLinearity[expr] //. {
  lf1 -> l1, sp[p1, p1] -> p^2, sp[l1, l1] -> l1^2,
  sp[l1, p1] -> l1*p*cos[l1, p], sp[p1, l1] -> l1*p*cos[l1, p],
  Sqrt[(a_)^2] :> a, ((a_)^2)^((n_)/2) :> a^n, cos[l1, p] :> cos1,
  Sqrt[(a_)^2]^(-1) :> a^(-1)};

(* vacuum multi-angle reduction for 3-/4-point vertices: sp/cos -> {l1, p, cosl1p_i} *)
SPParam[expr_] := UseLorentzLinearity[expr] //. {
  lf1 -> l1, sp[p, p] -> p^2, sp[l1, l1] -> l1^2,
  sp[l1, p1] -> p*l1*cos[l1, p1], sp[l1, p2] -> p*l1*cos[l1, p2],
  sp[l1, p3] -> p*l1*cos[l1, p3], sp[l1, p4] -> p*l1*cos[l1, p4],
  Sqrt[(a_)^2] :> a, ((a_)^2)^((n_)/2) :> a^n, ((l1_)^2)^(Rational[n_, 2]) :> l1^n,
  cos[l1, p1] :> cosl1p1, cos[l1, p2] :> cosl1p2, cos[l1, p3] :> cosl1p3, cos[l1, p4] :> cosl1p4};

(* Compile-time SU(N) sizes (baked, not runtime parameters). *)
SetNc[3];
SetNf[2];
$Assumptions = k > 0 && p > 0 && l1 > 0 && -1 < cos1 < 1 && -1 < cos2 < 1 && -1 < cos3 < 1;

(* ---- DiFfRG kernel parameter list (QCD.nb cell 6) ---------------------------------------------- *)
interpolatorType = "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>";
$ntDressings = {ZA3, ZAcbc, ZA4, ZAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, hPiL, hSigL, zPiL, zSigL};
kernelParameterList = Join[
  {<|"Name" -> "k", "Type" -> "double"|>, <|"Name" -> "etaPiL", "Type" -> "double"|>,
   <|"Name" -> "etaSigL", "Type" -> "double"|>, <|"Name" -> "d1V", "Type" -> "double"|>,
   <|"Name" -> "d2V", "Type" -> "double"|>, <|"Name" -> "rhoL", "Type" -> "double"|>},
  Table[<|"Name" -> ToString[d], "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>, {d, $ntDressings}]];

FSetRegisterSize[64];

(* Frame definitions (QCD.nb cell 7) live in the NUMERIC generator script — they use NumTracer's
   propFrame head, which is only defined once NumTracer.m is loaded (this setup is Get-ed first). *)
args0 = {l1, cos1, k};
args1 = {l1, cos1, p, k};

(* Codegen opt-in for the TESTS (as in qcd_setup.m): emit against the standalone shim so the numeric
   kernel and the copied FormTracer oracle share the DiFfRG namespace + complex/compute API. Each
   generator script applies it via SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts] after loading. *)
$ntShimOpts = {"RuntimeInclude" -> "shim.hpp", "ExtraIncludes" -> {},
  "KernelNamespace" -> "DiFfRG", "SupportNamespace" -> "DiFfRG",
  "DressingType" -> "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>"};
