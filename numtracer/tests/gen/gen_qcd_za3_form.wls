#!/usr/bin/env wolframscript
(* FormTracer (DiFfRG FORM-codegen) clean-cache TIMING for the QCD three-gluon vertex (ZA3),
   replicating QCD_vacuum_base/QCD.m's ZA3 block (2 angles, SP3 frame). Fresh empty TraceCache so
   FormTrace recomputes from scratch; TBCache stays warm. Times FORM trace + MakeKernel codegen.
   Run from numtracer/:   wolfram -script tests/gen/gen_qcd_za3_form.wls                         *)

If[Quiet[Needs["FunKit`"]] === $Failed, (* installed (Applications dir) or the local checkout *)
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "gen", "gen_qcd_form_setup.m"}]];

(* clean cold run: empty this kernel's trace cache so FormTrace reruns tform (the user's
   "EMPTY tracecache" requirement); TBCache is shared setup and left warm. *)
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_za3form"}];
If[DirectoryQ[$cacheDir] && !StringQ[Environment["NT_KEEP_TRACE_CACHE"]],
  DeleteDirectory[$cacheDir, DeleteContents -> True]]; (* NT_KEEP_TRACE_CACHE=1 reuses the cached tform run *)
FSetCacheDirectory[$cacheDir];

{tTrace, FlowA3} = AbsoluteTiming[
  Module[{fRG, proj, tr},
    fRG  = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3]}] // FTruncate // FPlot // FRoute // FPrint;
    proj = TBGetProjector["AAAClassTrans", 1, {i1, i2, i3} /. fRG["1-Loop"]["ExternalIndices"]] //
             TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3] & // Simplify;
    tr   = FTerm[proj] ** (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
    (FormTrace["ZA3", tr, {}, SP3FormRule] // dressingRules //
       TBProjectToSymmetricPoint[#, l1, p, p1, p2, p3] &) // SPParam // Simplify]];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..SPParam): ", tTrace, " s"];


(* WORKAROUND (2026-08-08): COEN's sub-kernel split emits `_cseN` references across sub-kernel
   scopes on the larger FORM kernels (undeclared-identifier C++, or a hard "cannot lift" abort).
   Disable the split — one flat kernel, same values, the same no-split fallback Cpp.m uses
   internally. Remove once the COEN lift handles renumbered `_cse` names across blocks. *)
FunKit`Private`$codeMaxKernelTerms = 10^9;
FunKit`Private`$kernelSplitLevel = 1;

tKernel = First @ AbsoluteTiming[
  MakeKernel[FlowA3,
    "Name" -> "ZA3", "Integrator" -> "Integrator_p2_4D_2ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> kernelParameterList, "KernelBody" -> SP3Defs,
    "IntegrationVariables" -> {"l1", "cos1", "cos2"},
    "Coordinates" -> {"LogarithmicCoordinates1D<double>"}, "CoordinateArguments" -> {"p"}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];
Print["[time] FormTracer generation TOTAL: ", tTrace + tKernel, " s"];

(* place the generated kernel header where the compare test includes it: flows/<Name>/kernel.hh
   (the pointwise integrand class — NOT the thin <Name>.hh integrator wrapper), with the DiFfRG
   framework includes swapped for the standalone shim. *)
$src = FileNameJoin[{$repo, "flows", "ZA3", "kernel.hh"}];
$out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "ZA3_kernel.hh"}];
If[FileExistsQ[$src],
  Module[{txt = Import[$src, "Text"]},
    txt = StringReplace[txt, RegularExpression["(?:#include \"DiFfRG/[^\"]*\"\n)+"] -> "#include \"shim.hpp\"\n"];
    Export[$out, txt, "Text"];
    Print["wrote oracle: ", $out]],
  Print["WARN: ", $src, " not found"]];
Print["DONE ZA3 FORM kernel."];
