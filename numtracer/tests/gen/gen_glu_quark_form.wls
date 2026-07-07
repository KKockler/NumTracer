#!/usr/bin/env wolframscript
(* FormTracer (FORM) oracle for the gluon <- quark loop flow (Glu_quark) — the independent
   cross-check for the NUMERIC matrix-product kernel, re-anchoring the compare test off the
   removed dense oracle. Builds the SAME hand-rolled quark-loop gluon self-energy that
   gen_glu_quark.wls builds (massless q-slash propagator, i*gamma^mu*T^a quark-gluon vertex,
   transverse/adjoint external projection), but traces it with FORM (FormTrace) instead of the
   numeric backend, then emits the kernel via DiFfRG's MakeKernel with a parameter list matching
   the numeric Glu_quark signature exactly:  kernel(l1, cos1, p, k)  — l1,cos1 the integration
   variables, p,k plain doubles, no dressings (RF/RFdot ride the Regulator template).  Baked
   SU(3) colour (Nc=3) and the closed quark flavour loop folded to Nf=2, matching the numeric
   SetNc[3]/SetNf[2].  Emits Glu_quark_form_kernel.hh into refshim/.  Run from numtracer/:
     wolfram -script tests/gen/gen_glu_quark_form.wls                                          *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];

(* Minimal truncation identical to gen_glu_quark.wls: only the quark loop closes the gluon
   two-point, so FTakeDerivatives[.., {A, A}] yields exactly the quark-loop diagram. *)
fields = <|"Commuting" -> {A[p, {v, c}]}, "Grassmann" -> {{qb[p, {d, A, F}], q[p, {d, A, F}]}}|>;
truncation = <|GammaN -> {{A, A}, {A, qb, q}, {qb, q}}, Propagator -> {{A, A}, {qb, q}},
  R -> {{qb, q}}, Rdot -> {{qb, q}}, Field -> {{}}|>;
FSetGlobalSetup[<|"FieldSpace" -> fields, "Truncation" -> truncation,
  "FeynmanRules" -> <|GammaN -> {}, Propagator -> {}, R -> {}, Rdot -> {}|>|>];

(* clean cold run: empty this kernel's trace cache so FormTrace reruns tform; TBCache stays warm. *)
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_glu_quark_form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

fRG = FTakeDerivatives[WetterichEquation, {A[i1], A[i2]}] // FTruncate // FRoute;
expr = fRG["1-Loop"]["Expression"];

(* SAME hand-rolled diagrammatic rules as gen_glu_quark.wls. QG vertex: i gamma^mu (Dirac+Lorentz)
   x T^a (fundamental colour) x delta (trivial flavour). Quark prop / Rdot: massless q-slash x
   delta_color x delta_flav, with the propagator/Rdot dressings kept ABSTRACT (dressing[...] symbols)
   so FORM traces them as scalar coefficients — substituted (massless simplification) after the
   trace, mirroring the Zq FORM oracle's post-trace dressingRules. *)
myRules = {
  GammaN[{A, qb, q}, {{_, {lor_, adjc_}}, {_, {dqb_, cqb_, fqb_}}, {_, {dq_, cq_, fq_}}}] :>
    I gamma[lor, dqb, dq] TCol[adjc, cqb, cq] deltaFundFlav[fqb, fq],
  Propagator[{q, qb} | {qb, q}, {{m1_, {d1_, c1_, f1_}}, {m2_, {d2_, c2_, f2_}}}] :>
    With[{mu = Unique["lm"]},
      deltaFundCol[c1, c2] deltaFundFlav[f1, f2] *
      (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[GammaN, {qb, q}, 1, {m1, m2}]],
  Rdot[{qb, q}, {{m1_, {d1_, c1_, f1_}}, {m2_, {d2_, c2_, f2_}}}] :>
    With[{mu = Unique["lm"]},
      deltaFundCol[c1, c2] deltaFundFlav[f1, f2] *
      (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[Rdot, {qb, q}, 1, {m1, m2}]]};

(* project the external gluon (i1,i2): transverse projector / (d-1) on Lorentz, adjoint delta /
   (Nc^2-1) on colour;  1/((4-1)(3^2-1)) = 1/24. *)
trace = FTerm[1/24, ntTransProjTok[p1, v1, v2] deltaAdjCol[c1, c2]] ** (expr /. myRules) /.
  ntTransProjTok -> transProj;

(* FORM trace: Dirac gamma chain, fundamental colour tr(T^a T^b) (Nc=3 baked), trivial closed
   flavour loop (deltaFundFlav -> Nf).  Keeps the abstract dressing[...] symbols. *)
{tTrace, traced} = AbsoluteTiming[FormTrace["Glu_quark_form", trace]];
Print["[time] FORM trace: ", tTrace, " s"];

(* massless simplified dressings, SAME as gen_glu_quark.wls: 1/(q^2 + RF), RFdot. *)
dressed = traced //. {
  dressing[GammaN, {qb, q}, 1, {m1_, m2_}] :> 1/(sp[m2, m2] + RF[k^2, sp[m2, m2]]),
  dressing[Rdot, {qb, q}, 1, {m1_, m2_}] :> RFdot[k^2, sp[m2, m2]]};

(* vacuum 1-angle scalar reduction, matching the numeric frame propFrame[p,l1,cos1,p1,l1]. *)
PropParam[e_] := UseLorentzLinearity[e] //. {
  lf1 -> l1,
  sp[p1, p1] -> p^2, sp[l1, l1] -> l1^2,
  sp[l1, p1] -> l1 p cos1, sp[p1, l1] -> l1 p cos1,
  Sqrt[a_^2] :> a, (a_^2)^(n_/2) :> a^n};

(* fold the closed quark flavour loop to Nf = 2 (numeric bakes SetNf[2]); colour Nc = 3 baked. *)
flow = (Total @ Flatten @ {dressed}) // PropParam;
flow = Simplify[flow /. {Global`Nf -> 2, Global`Nc -> 3}];

tKernel = First @ AbsoluteTiming[
  MakeKernel[flow,
    "Name" -> "Glu_quark_form", "Integrator" -> "Integrator_p2_1ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> {<|"Name" -> "p", "Type" -> "double"|>, <|"Name" -> "k", "Type" -> "double"|>},
    "IntegrationVariables" -> {"l1", "cos1"},
    "Coordinates" -> {}, "CoordinateArguments" -> {}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];

(* place the generated kernel where the compare test includes it: swap DiFfRG physics includes for
   the standalone shim, matching the other refshim/*_form_kernel.hh oracles. *)
$src = FileNameJoin[{$repo, "flows", "Glu_quark_form", "kernel.hh"}];
$out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "Glu_quark_form_kernel.hh"}];
If[FileExistsQ[$src],
  Module[{txt},
    txt = Import[$src, "Text"];
    txt = StringReplace[txt,
      "#include \"DiFfRG/physics/interpolation.hh\"\n#include \"DiFfRG/physics/physics.hh\"" ->
        "#include \"shim.hpp\""];
    txt = StringReplace[txt, "#include \"DiFfRG/physics/physics.hh\"" -> "#include \"shim.hpp\""];
    Export[$out, txt, "Text"];
    Print["wrote oracle: ", $out]],
  Print["WARN: ", $src, " not found"]];
Print["DONE Glu_quark FORM kernel."];
