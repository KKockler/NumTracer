#!/usr/bin/env wolframscript
(* FormTracer (FORM) oracle for the pure-YM ZA4 four-gluon vertex flow — the independent
   cross-check for gen_ym_za4_numeric.wls. Traces the SAME loop (A^4 derivative of the Wetterich
   equation, projected with TBGetProjector["AAAAClassTrans", 1] at the symmetric point) with FORM,
   and emits the kernel via DiFfRG's MakeKernel with a parameter list matching the numeric YM ZA4
   kernel's EXACT signature so the validation is a pointwise compare.
   Emits refshim/YM_ZA4_kernel.hh. Regenerate this TOGETHER with gen_ym_za4_numeric.wls and commit
   both in one commit — the pre-existing refshim/ZA4_kernel.hh vs gen/ZA4_num_kernel.hh pair drifted
   to different vintages, which silently degraded that comparison. *)
If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "ym_setup.m"}]];

(* clean cold run: a cached trace is served silently and fakes agreement *)
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_ym_za4_form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

interpolatorType = "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>";
za4Dressings = {ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA};
za4ParameterList = Join[
  {<|"Name" -> "p", "Type" -> "double"|>, <|"Name" -> "k", "Type" -> "double"|>},
  Table[<|"Name" -> ToString[d], "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>,
    {d, za4Dressings}]];

(* Full permutation group of the identical external legs; the ClassTrans/transAAqbq
   projectors are covariant under it. Symmetries must be stated by hand since FunKit
   0bacdd1 -- without this the derivative yields every leg ordering as its own diagram
   (measured on za4_147: 6183 diagrams / 164 s, vs 294 / 7.6 s with it). *)
symLegs = FMakeSymmetryList[{A[i1], A[i2], A[i3], A[i4]}];
fRGA4 = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3], A[i4]}, "Symmetries" -> symLegs] // FTruncate // FPlot // FRoute // FPrint;
projectorA4 = FTerm[TBGetProjector["AAAAClassTrans", 1, {i1, i2, i3, i4} /. fRGA4["1-Loop"]["ExternalIndices"]]];
traceExprA4 = projectorA4 ** (fRGA4["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
{tTrace, FlowA4} = AbsoluteTiming[
  FormTrace["YM_ZA4_form", traceExprA4, {}, SP4FormRule] // dressingRules //
    TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3, p4] & // SPParam // Simplify];
Print["[time] FORM trace: ", tTrace, " s"];

tKernel = First @ AbsoluteTiming[
  MakeKernel[FlowA4,
    "Name" -> "YM_ZA4_form", "Integrator" -> "Integrator_p2_4D_3ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> za4ParameterList, "KernelBody" -> SP4Defs,
    "IntegrationVariables" -> {"l1", "cos1", "cos2", "phi"},
    "Coordinates" -> {}, "CoordinateArguments" -> {}]];
Print["[time] FormTracer MakeKernel: ", tKernel, " s"];

$src = FileNameJoin[{$repo, "flows", "YM_ZA4_form", "kernel.hh"}];
$out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "YM_ZA4_kernel.hh"}];
If[FileExistsQ[$src],
  Module[{txt = Import[$src, "Text"]},
    txt = StringReplace[txt,
      "#include \"DiFfRG/physics/interpolation.hh\"\n#include \"DiFfRG/physics/physics.hh\"" ->
        "#include \"shim.hpp\""];
    Export[$out, txt, "Text"];
    Print["wrote oracle: ", $out]],
  Print["WARN: ", $src, " not found"]];
Print["DONE YM ZA4 FORM kernel."];
