#!/usr/bin/env wolframscript
(* Pure-YM ZA4 (four-gluon vertex) via the NUMERIC matrix-product backend.
   This is the flow measured broken in NUMTRACER-ZA4-vertex-bug.md: its ZA3^2.ZA4 triangle and
   ZA4^2 bubble disagree with FORM at ~0.55 scale-relative while its ZA3^4 box agrees to 1.3e-14.
   The pre-existing ZA4 coverage (gen_qcd_za4_numeric.wls) is the QCD *quark box*, which shares
   none of those topologies — hence this generator. Oracle: gen_ym_za4_form.wls. *)
If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "ym_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];

frame = NumTracer`sp4Frame[Global`p, Global`l1, Global`cos1, Global`cos2, Global`phi,
          Global`p1, Global`p2, Global`p3, Global`p4, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`phi, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, dtZc, Zc, dtZA, ZA};
spAngles4 = {
  Global`cosl1p1 -> Sqrt[1 - Global`cos1^2] Global`cos2,
  Global`cosl1p2 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 - 2 Sqrt[2 - 2 Global`cos2^2] Cos[Global`phi]),
  Global`cosl1p3 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 + Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] - Sqrt[3] Sin[Global`phi])),
  Global`cosl1p4 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 + Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] + Sqrt[3] Sin[Global`phi]))};

fRG  = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3], A[i4]}] // FTruncate // FPlot // FRoute // FPrint;
proj = FTerm[TBGetProjector["AAAAClassTrans", 1, {i1, i2, i3, i4} /. fRG["1-Loop"]["ExternalIndices"]]];
traceA4 = proj ** (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
net = FromFunKit[ SPParam[traceA4 // dressingRules] ];
Print["YM ZA4 net top terms = ", If[Head[net] === Plus, Length[net], 1]];

{tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
Print["[time] NumTrace (analyse): ", tA, " s   diagrams = ", Length[ntk[[1]]["Diagrams"]]];

tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_ym_za4_num.cpp"}], FileNameJoin[{$gen, "YM_ZA4_num_kernel.hh"}],
  FileNameJoin[{$gen, "YM_ZA4_num_kernels.hh"}],
  "Name" -> "YM_ZA4_num_kernel", "Namespace" -> "ym_za4_num", "Dressings" -> dress,
  "AngleDefs" -> spAngles4] ];
Print["[time] Numeric backend: ", tN, " s"];
Print["YM ZA4 numeric kernel generated."];
