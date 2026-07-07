#!/usr/bin/env wolframscript
(* FormTracer (FORM) oracle for the pure-YM ZAcbc ghost-gluon vertex flow — the independent
   cross-check for the NUMERIC matrix-product kernel (re-anchoring the ZAcbc test off the removed
   Dense oracle). Traces the SAME ghost-gluon two-point loop as gen_ym_zacbc_numeric.wls
   (A-cb-c derivative of the Wetterich equation, projected with TBGetProjector["Acbc", 1] at the
   symmetric point) with FORM instead of the numeric backend, and emits the kernel via DiFfRG's
   MakeKernel with a parameter list matching the numeric ZAcbc kernel's EXACT signature
   (baked SU(3), same ordered dressings, p and k as plain doubles) so the validation is a pointwise
   compare. Emits ZAcbc_form_kernel.hh into refshim/. Run from numtracer/:
     wolfram -script tests/gen/gen_zacbc_form.wls                                                 *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "ym_setup.m"}]];  (* bakes SU(3); dressingRules / SPParam / SP3FormRule / SP3Defs *)

(* clean cold run: empty this kernel's trace cache so FormTrace reruns tform; TBCache stays warm. *)
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_zacbc_form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

(* DiFfRG kernel parameter list matching the NUMERIC ZAcbc signature (gen_ym_zacbc_numeric.wls):
   kernel(l1, cos1, cos2, p, k, ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA). l1, cos1, cos2 are integration
   variables; p, k are plain doubles (so Coordinates/CoordinateArguments are empty and p is a
   Parameter, unlike Yang-Mills.m which makes p a coordinate); the seven dressings are shim
   SplineInterpolator1D references, in the numeric kernel's order. *)
interpolatorType = "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>";
zacbcDressings = {ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA};
zacbcParameterList = Join[
  {<|"Name" -> "p", "Type" -> "double"|>, <|"Name" -> "k", "Type" -> "double"|>},
  Table[<|"Name" -> ToString[d], "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>, {d, zacbcDressings}]];

(* ---- ZAcbc ghost-gluon vertex loop (same construction as Yang-Mills.m / gen_ym_zacbc_numeric.wls) ---- *)
fRGAcbc = FTakeDerivatives[WetterichEquation, {A[i1], cb[i2], c[i3]}] // FTruncate // FSimplify // FRoute;
projectorAcbc = FTerm[TBGetProjector["Acbc", 1, {i1, i2, i3} /. fRGAcbc["1-Loop"]["ExternalIndices"]]];
traceExprAcbc = projectorAcbc ** (fRGAcbc["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
{tTrace, FlowAcbc} = AbsoluteTiming[
  FormTrace["ZAcbc_form", traceExprAcbc, {}, SP3FormRule] // dressingRules //
    TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3] & // SPParam // Simplify];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..SPParam): ", tTrace, " s"];

tKernel = First @ AbsoluteTiming[
  MakeKernel[FlowAcbc,
    "Name" -> "ZAcbc_form", "Integrator" -> "Integrator_p2_4D_2ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> zacbcParameterList, "KernelBody" -> SP3Defs,
    "IntegrationVariables" -> {"l1", "cos1", "cos2"},
    "Coordinates" -> {}, "CoordinateArguments" -> {}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];

(* place the generated kernel header where the compare test includes it: swap DiFfRG physics
   includes for the standalone shim, matching the other refshim/*_form_kernel.hh oracles. *)
$src = FileNameJoin[{$repo, "flows", "ZAcbc_form", "kernel.hh"}];
$out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "ZAcbc_form_kernel.hh"}];
If[FileExistsQ[$src],
  Module[{txt},
    txt = Import[$src, "Text"];
    txt = StringReplace[txt,
      "#include \"DiFfRG/physics/interpolation.hh\"\n#include \"DiFfRG/physics/physics.hh\"" ->
        "#include \"shim.hpp\""];
    Export[$out, txt, "Text"];
    Print["wrote oracle: ", $out]],
  Print["WARN: ", $src, " not found"]];
Print["DONE ZAcbc FORM kernel."];
