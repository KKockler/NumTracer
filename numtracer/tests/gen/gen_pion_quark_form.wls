#!/usr/bin/env wolframscript
(* FormTracer (FORM) oracle for the Phase-2 pion<-quark-loop->pion Yukawa flow (Pion_quark) — the
   independent cross-check for the NUMERIC matrix-product kernel, re-anchoring the test off the
   removed Dense oracle. Builds the SAME one-loop pion self-energy that gen_pion_quark.wls builds
   (closed massless quark line, two i*gamma5*(tau*pi) Yukawa vertices, a fermion regulator R_F ~
   q-slash, baked closed-form fermion dressings G = 1/(q^2 + R_F) and R_Fdot), but traces it with
   FORM (colour SU(3) tr(1)=Nc, flavour SU(2) tr(tau^a tau^b), and the Dirac trace) instead of the
   numeric backend, then emits the kernel via DiFfRG's MakeKernel with a Parameters list matching
   the numeric Pion_quark kernel's EXACT signature kernel(l1, cos1, p, k) — baked SU(N), no dressing
   parameters. Emits Pion_quark_form_kernel.hh into refshim/. Run from numtracer/:
     wolfram -script tests/gen/gen_pion_quark_form.wls                                              *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];

(* clean cold run: empty this kernel's trace cache so FormTrace reruns tform; TBCache stays warm. *)
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_pion_quark_form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

(* Compile-time SU(N): colour SU(3), flavour SU(2). Baked here so the FORM kernel and the numeric
   kernel share one signature. *)
SetNc[3];
SetNf[2];
$Assumptions = k > 0 && p > 0 && l1 > 0 && -1 < cos1 < 1;

(* ---- field space + truncation (mirror of gen_pion_quark.wls) ----------------------------------- *)
fields = <|"Commuting" -> {Pi[p, {i}]}, "Grassmann" -> {{qb[p, {d, A, F}], q[p, {d, A, F}]}}|>;
truncation = <|GammaN -> {{Pi, Pi}, {Pi, qb, q}, {qb, q}}, Propagator -> {{Pi, Pi}, {qb, q}},
  R -> {{qb, q}}, Rdot -> {{qb, q}}, Field -> {{}}|>;
FSetGlobalSetup[<|"FieldSpace" -> fields, "Truncation" -> truncation,
  "FeynmanRules" -> <|GammaN -> {}, Propagator -> {}, R -> {}, Rdot -> {}|>|>];

fRG = FTakeDerivatives[WetterichEquation, {Pi[i1], Pi[i2]}] // FTruncate // FRoute;
expr = fRG["1-Loop"]["Expression"];

(* ---- hand-rolled rules, FormTracer-native heads. Yukawa: i gamma5 (Dirac) x tau^a (flavour
        generator TFlav) x delta (colour). Quark prop / Rdot: massless q-slash numerator x
        deltaFundCol (colour) x deltaFundFlav (in-trace fundamental flavour). ---- *)
myRules = {
  GammaN[{Pi, qb, q}, {{_, {iso_}}, {_, {dqb_, cqb_, fqb_}}, {_, {dq_, cq_, fq_}}}] :>
    I gamma5[dqb, dq] TFlav[iso, fqb, fq] deltaFundCol[cqb, cq],
  Propagator[{q, qb} | {qb, q}, {{m1_, {d1_, c1_, f1_}}, {m2_, {d2_, c2_, f2_}}}] :>
    With[{mu = Unique["lm"]},
      deltaFundCol[c1, c2] deltaFundFlav[f1, f2] *
      (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[GammaN, {qb, q}, 1, {m1, m2}]],
  Rdot[{qb, q}, {{m1_, {d1_, c1_, f1_}}, {m2_, {d2_, c2_, f2_}}}] :>
    With[{mu = Unique["lm"]},
      deltaFundCol[c1, c2] deltaFundFlav[f1, f2] *
      (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[Rdot, {qb, q}, 1, {m1, m2}]]};

(* project the external isospin (i1,i2) onto the scalar self-energy *)
trace = FTerm[1/3, deltaAdjFlav[Global`i1, Global`i2]] ** (expr /. myRules);

(* baked closed-form fermion dressings: G(q^2) = 1/(q^2 + RF(k^2,q^2)) ; regulator-dot RFdot. *)
dressed = trace //. {
  dressing[GammaN, {qb, q}, 1, {m1_, m2_}] :> 1/(sp[m2, m2] + RF[k^2, sp[m2, m2]]),
  dressing[Rdot, {qb, q}, 1, {m1_, m2_}] :> RFdot[k^2, sp[m2, m2]]};
dressed = UseLorentzLinearity[dressed] //. {lf1 -> l1};

(* ---- FORM trace + scalar (1-angle) reduction -------------------------------------------------- *)
(* 2-point momentum conservation p2 -> -p1, then external p1 -> magnitude p, loop l1, angle cos1. *)
PropParam[e_] := UseLorentzLinearity[e] //. {
  lf1 -> l1, sp[p1, p1] -> p^2, sp[l1, l1] -> l1^2,
  sp[l1, p1] -> l1*p*cos[l1, p], sp[p1, l1] -> l1*p*cos[l1, p],
  Sqrt[(a_)^2] :> a, ((a_)^2)^((n_)/2) :> a^n, cos[l1, p] :> cos1,
  Sqrt[(a_)^2]^(-1) :> a^(-1)};

{tTrace, Flow} = AbsoluteTiming[
  PropParam[(FormTrace["Pion_quark_form", dressed] // Simplify) //. Global`p2 -> -Global`p1]];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..PropParam): ", tTrace, " s"];
Print["Flow head = ", Head[Flow], "  (nonzero: ", Simplify[Flow] =!= 0, ")"];

(* DiFfRG kernel parameter list matching the NUMERIC Pion_quark signature (gen_pion_quark.wls):
   kernel(l1, cos1, p, k). l1, cos1 are integration variables; p, k are plain doubles; NO dressings. *)
pionParameterList = {<|"Name" -> "p", "Type" -> "double"|>, <|"Name" -> "k", "Type" -> "double"|>};

tKernel = First @ AbsoluteTiming[
  MakeKernel[Flow,
    "Name" -> "Pion_quark_form", "Integrator" -> "Integrator_p2_1ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> pionParameterList,
    "IntegrationVariables" -> {"l1", "cos1"},
    "Coordinates" -> {}, "CoordinateArguments" -> {}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];

(* place the generated kernel header where the compare test includes it: swap DiFfRG physics
   includes for the standalone shim, matching the other refshim/*_form_kernel.hh oracles. *)
$src = FileNameJoin[{$repo, "flows", "Pion_quark_form", "kernel.hh"}];
$out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "Pion_quark_form_kernel.hh"}];
If[FileExistsQ[$src],
  Module[{txt},
    txt = Import[$src, "Text"];
    txt = StringReplace[txt,
      "#include \"DiFfRG/physics/interpolation.hh\"\n#include \"DiFfRG/physics/physics.hh\"" ->
        "#include \"shim.hpp\""];
    Export[$out, txt, "Text"];
    Print["wrote oracle: ", $out]],
  Print["WARN: ", $src, " not found"]];
Print["DONE Pion_quark FORM kernel."];
