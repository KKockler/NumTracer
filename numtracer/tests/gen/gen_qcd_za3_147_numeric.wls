#!/usr/bin/env wolframscript
(* QCD ZA3 (three-gluon vertex, with quark triangle) using the FULL quark-gluon vertex {1,4,7}
   tensor basis, via the NUMERIC matrix-product backend — the dense `1/4/7` headline target.
   Sibling of gen_qcd_za3_numeric.wls (which is struct-1 only); loads qcd147_setup.m and emits
   distinct kernel names so the struct-1 fixture is untouched. Run from numtracer/:
     wolfram -script tests/gen/gen_qcd_za3_147_numeric.wls                                        *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd147_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];
(* opt IN to vertex collection (default OFF since it OOMs high-P flows): za3_147 is the small-P flow
   where collecting the full 1/4/7 vertex wins (8.9 s vs 13.3 s) and it exercises the ntDiracSlot path *)
NumTracer`Private`$ntVertexCollect = True;
(* the full 1/4/7 vertex basis emits many more (per-branch) numeric entries than struct-1; the
   codegen expression tree exceeds Mathematica's conservative default recursion depth (1024). *)
$RecursionLimit = 1000000;

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
          Global`p1, Global`p2, Global`p3, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4, ZAqbq7, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};

(* AAAClassTrans is Bose-symmetric in all three gluon legs at generic momenta, so the FULL
   permutation group may be used — the projector below is covariant under every element. Since
   FunKit 0bacdd1 symmetries must be stated by hand; without this the derivative produces every
   leg ordering as its own diagram. Same choice as production QCD.wl. *)
symA3 = FMakeSymmetryList[{A[i1], A[i2], A[i3]}];
fRGA3 = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3]}, "Symmetries" -> symA3] // FTruncate // FPlot // FRoute // FPrint;
projA3 = FTerm[TBGetProjector["AAAClassTrans", 1, {i1, i2, i3} /. fRGA3["1-Loop"]["ExternalIndices"]]];
traceA3 = projA3 ** (fRGA3["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
spAngles = {
  Global`cosl1p1 -> Global`cos1,
  Global`cosl1p2 -> (-Global`cos1 + Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2,
  Global`cosl1p3 -> (-Global`cos1 - Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2};
flow = SPParam[traceA3 // dressingRules];
net  = FromFunKit[flow];
Print["ZA3 147 net top terms = ", If[Head[net] === Plus, Length[net], 1]];

{tA, ntk} = AbsoluteTiming[ NumTrace[net, "Frame" -> frame, "Args" -> args] ];
Print["[time] NumTrace (analyse): ", tA, " s   diagrams = ", Length[ntk[[1]]["Diagrams"]]];

tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_qcd_za3_147_num.cpp"}],
  FileNameJoin[{$gen, "ZA3_147_num_kernel.hh"}],
  FileNameJoin[{$gen, "ZA3_147_num_kernels.hh"}],
  "Name" -> "ZA3_147_num_kernel", "Namespace" -> "za3_147_num",
  "Dressings" -> dress, "AngleDefs" -> spAngles] ];
Print["[time] Numeric backend (contract+compile+run): ", tN, " s"];
Print["ZA3 147 numeric kernels generated."];
