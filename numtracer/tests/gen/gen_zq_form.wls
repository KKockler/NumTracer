#!/usr/bin/env wolframscript
(* FormTracer (FORM) oracle for the QCD quark wave-function Zq — the independent cross-check for the
   NUMERIC matrix-product kernel. Traces the SAME quark self-energy loop (qbq derivative of the
   Wetterich equation, projected with TBGetProjector["qbq",1]) with FORM instead of the numeric
   backend, and emits the kernel via DiFfRG's MakeKernel with a parameter list matching the numeric
   Zq kernel (baked SU(3)/SU(2), same ordered dressings) so the validation is a pointwise compare.
   Emits Zq_form_kernel.hh into refshim/. Run from numtracer/:
     wolfram -script tests/gen/gen_zq_form.wls                                                      *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_setup.m"}]];  (* bakes SU(3)/SU(2) *)

(* clean cold run: empty this kernel's trace cache so FormTrace reruns tform; TBCache stays warm. *)
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_zq_form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

(* DiFfRG kernel parameter list matching the NUMERIC Zq signature (gen_zq_numeric.wls):
   kernel(l1, cos1, p, k, ZA, Zq, Mq, ZAqbq1, dtZA, dtZq). l1, cos1 are integration variables;
   p, k are plain doubles; the six dressings are shim SplineInterpolator1D references. *)
interpolatorType = "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>";
zqDressings = {ZA, Zq, Mq, ZAqbq1, dtZA, dtZq};
zqParameterList = Join[
  {<|"Name" -> "p", "Type" -> "double"|>, <|"Name" -> "k", "Type" -> "double"|>},
  Table[<|"Name" -> ToString[d], "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>, {d, zqDressings}]];

(* ---- Zq quark self-energy loop (same construction as gen_zq_numeric.wls / gen_qcd_zq_form.wls) --- *)
fRGqbq = FTakeDerivatives[WetterichEquation, {qb[i1], q[i2]}] // FTruncate // FRoute;
tr = FTerm[-TBGetProjector["qbq", 1, {i2, i1} /. fRGqbq["1-Loop"]["ExternalIndices"]]] **
       (fRGqbq["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
{tTrace, Flowqbq} = AbsoluteTiming[(FormTrace["Zq_form", tr] // dressingRules) // PropParam // Simplify];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..PropParam): ", tTrace, " s"];

tKernel = First @ AbsoluteTiming[
  MakeKernel[Flowqbq,
    "Name" -> "Zq_form", "Integrator" -> "Integrator_p2_1ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> zqParameterList,
    "IntegrationVariables" -> {"l1", "cos1"},
    "Coordinates" -> {}, "CoordinateArguments" -> {}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];

(* place the generated kernel header where the compare test includes it: swap DiFfRG physics
   includes for the standalone shim, matching the other refshim/*_form_kernel.hh oracles. *)
$src = FileNameJoin[{$repo, "flows", "Zq_form", "kernel.hh"}];
$out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "Zq_form_kernel.hh"}];
If[FileExistsQ[$src],
  Module[{txt},
    txt = Import[$src, "Text"];
    txt = StringReplace[txt,
      "#include \"DiFfRG/physics/interpolation.hh\"\n#include \"DiFfRG/physics/physics.hh\"" ->
        "#include \"shim.hpp\""];
    Export[$out, txt, "Text"];
    Print["wrote oracle: ", $out]],
  Print["WARN: ", $src, " not found"]];
Print["DONE Zq FORM kernel."];
