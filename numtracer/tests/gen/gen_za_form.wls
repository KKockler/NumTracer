#!/usr/bin/env wolframscript
(* FormTracer (FORM) oracle for the QCD gluon self-energy ZA — the independent cross-check for the
   NUMERIC matrix-product kernel. Traces the SAME gluon self-energy loop (AA derivative of the
   Wetterich equation, projected with TBGetProjector["AA",1] and the 1/sp[p1,p1] normalisation) with
   FORM instead of the numeric backend, and emits the kernel via DiFfRG's MakeKernel with a parameter
   list matching the numeric ZA kernel (baked SU(3)/SU(2), same ordered dressings) so the validation is
   a pointwise / cos1-integrated compare. Emits ZA_form_kernel.hh into refshim/. Run from numtracer/:
     wolfram -script tests/gen/gen_za_form.wls                                                       *)

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
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_za_form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

(* DiFfRG kernel parameter list matching the NUMERIC ZA signature (gen_qcd_za_numeric.wls):
   kernel(l1, cos1, p, k, ZA3, ZAcbc, ZA4, ZAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq). l1, cos1 are
   integration variables; p, k are plain doubles; the eleven dressings are shim SplineInterpolator1D
   references (Nf baked, so NO Nf parameter). *)
interpolatorType = "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>";
zaDressings = {ZA3, ZAcbc, ZA4, ZAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};
zaParameterList = Join[
  {<|"Name" -> "p", "Type" -> "double"|>, <|"Name" -> "k", "Type" -> "double"|>},
  Table[<|"Name" -> ToString[d], "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>, {d, zaDressings}]];

(* ---- ZA gluon self-energy loop (same construction as gen_qcd_za_numeric.wls) ---- *)
fRGAA = FTakeDerivatives[WetterichEquation, {A[i1], A[i2]}] // FTruncate // FRoute;
tr = FTerm[1/sp[p1, p1], TBGetProjector["AA", 1, {i1, i2} /. fRGAA["1-Loop"]["ExternalIndices"]]] **
       (fRGAA["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]) //. p2 -> -p1;
{tTrace, FlowAA} = AbsoluteTiming[(FormTrace["ZA_form", tr] // dressingRules) // PropParam // Simplify];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..PropParam): ", tTrace, " s"];

tKernel = First @ AbsoluteTiming[
  MakeKernel[FlowAA,
    "Name" -> "ZA_form", "Integrator" -> "Integrator_p2_1ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> zaParameterList,
    "IntegrationVariables" -> {"l1", "cos1"},
    "Coordinates" -> {}, "CoordinateArguments" -> {}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];

(* place the generated kernel header where the compare test includes it: swap DiFfRG physics
   includes for the standalone shim, matching the other refshim/*_form_kernel.hh oracles. *)
$src = FileNameJoin[{$repo, "flows", "ZA_form", "kernel.hh"}];
$out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "ZA_form_kernel.hh"}];
If[FileExistsQ[$src],
  Module[{txt},
    txt = Import[$src, "Text"];
    txt = StringReplace[txt,
      "#include \"DiFfRG/physics/interpolation.hh\"\n#include \"DiFfRG/physics/physics.hh\"" ->
        "#include \"shim.hpp\""];
    Export[$out, txt, "Text"];
    Print["wrote oracle: ", $out]],
  Print["WARN: ", $src, " not found"]];
Print["DONE ZA FORM kernel."];
