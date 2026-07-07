#!/usr/bin/env wolframscript
(* FormTracer (DiFfRG FORM-codegen) oracle for the QCD pion anomalous dimension etaPiL — the
   independent cross-check for the NUMERIC matrix-product kernel. Traces the SAME pion two-point
   loop (QCD.nb cells 14-15) with FORM instead of the numeric backend, and emits the kernel via
   DiFfRG's MakeKernel. Both backends bake the SU(N) sizes (colour SU(3), flavour SU(2)) so the FORM
   kernel and the numeric kernel share one signature, making the validation a pointwise compare.
   Emits EtaPiL_form_kernel.hh into refshim/. Run from numtracer/:
     wolfram -script tests/gen/gen_qcd_etapil_form.wls                                              *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_meson_setup.m"}]];

(* clean cold run: empty this kernel's trace cache so FormTrace reruns tform; TBCache stays warm. *)
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_etapil_form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

(* ---- pion two-point loop + wave-function extraction (QCD.nb cells 14-15) ----------------------- *)
fRGPiPi = FTakeDerivatives[WetterichEquation, {ΠL[i1], ΠL[i2]}] // FTruncate // FRoute;
traceΠL = FTerm[deltaAdjFlav[f1, f2], fRGPiPi["1-Loop"]["Expression"] /. diagRules] //. p2 -> -p1;
traceΠL = (Times @@ #1 &) /@ List @@ traceΠL;
traceΠL = Simplify[UseLorentzLinearity[(traceΠL - (traceΠL //. p1 -> 0))/sp[p1, p1]]];
traceΠL = FEx @@ FTerm /@ traceΠL;

{tTrace, FlowPiL} = AbsoluteTiming[
  (FormTrace["etaPiL", FTerm[-1] ** traceΠL] // dressingRules // PropParam // Simplify) /. p -> 10^(-3)];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..PropParam): ", tTrace, " s"];

tKernel = First @ AbsoluteTiming[
  MakeKernel[FlowPiL,
    "Name" -> "EtaPiL_form", "Integrator" -> "Integrator_p2_1ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> kernelParameterList,
    "IntegrationVariables" -> {"l1", "cos1"},
    "Coordinates" -> {}, "CoordinateArguments" -> {}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];

(* place the generated kernel header where the compare test includes it: take the FORM kernel class
   (flows/EtaPiL_form/kernel.hh) and swap its DiFfRG physics includes for the standalone shim, so the
   oracle compiles dependency-free exactly like the other refshim/*_kernel.hh FormTracer oracles. *)
$src = FileNameJoin[{$repo, "flows", "EtaPiL_form", "kernel.hh"}];
$out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "EtaPiL_form_kernel.hh"}];
If[FileExistsQ[$src],
  Module[{txt},
    txt = Import[$src, "Text"];
    txt = StringReplace[txt,
      "#include \"DiFfRG/physics/interpolation.hh\"\n#include \"DiFfRG/physics/physics.hh\"" ->
        "#include \"shim.hpp\""];
    Export[$out, txt, "Text"];
    Print["wrote oracle: ", $out]],
  Print["WARN: ", $src, " not found"]];
Print["DONE etaPiL FORM kernel."];
