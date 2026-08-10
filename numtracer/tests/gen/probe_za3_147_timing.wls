#!/usr/bin/env wolframscript
(* Stage 4: generation timing, collection ON vs OFF, on za3_147 (full AqbqDirect147 vertex).
   Writes to throwaway paths so the committed collected kernel is untouched. *)
If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
$tmp  = "/tmp/za3timing";
Quiet@CreateDirectory[$tmp];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd147_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];
$RecursionLimit = 1000000;
frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2, Global`p1, Global`p2, Global`p3, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4, ZAqbq7, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};
spAngles = {Global`cosl1p1 -> Global`cos1,
  Global`cosl1p2 -> (-Global`cos1 + Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2,
  Global`cosl1p3 -> (-Global`cos1 - Sqrt[3 - 3 Global`cos1^2] Global`cos2)/2};
(* Full permutation group of the identical external legs; the ClassTrans/transAAqbq
   projectors are covariant under it. Symmetries must be stated by hand since FunKit
   0bacdd1 -- without this the derivative yields every leg ordering as its own diagram
   (measured on za4_147: 6183 diagrams / 164 s, vs 294 / 7.6 s with it). *)
symLegs = FMakeSymmetryList[{A[i1], A[i2], A[i3]}];
fRGA3 = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3]}, "Symmetries" -> symLegs] // FTruncate // FRoute;
projA3 = FTerm[TBGetProjector["AAAClassTrans", 1, {i1, i2, i3} /. fRGA3["1-Loop"]["ExternalIndices"]]];
traceA3 = projA3 ** (fRGA3["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
flow = SPParam[traceA3 // dressingRules];

timeGen[coll_, tag_, ns_] := Module[{net, tA, ntk, tG},
  net = FromFunKit[flow, "DressingCollection" -> coll];
  {tA, ntk} = AbsoluteTiming[NumTrace[net, "Frame" -> frame, "Args" -> args, "DressingCollection" -> coll]];
  tG = First @ AbsoluteTiming[ MakeNTKernel[ntk,
     FileNameJoin[{$tmp, "gen_" <> ns <> ".cpp"}],
     FileNameJoin[{$tmp, ns <> "_kernel.hh"}],
     FileNameJoin[{$tmp, ns <> "_kernels.hh"}],
     "Name" -> ns <> "_kernel", "Namespace" -> ns, "Dressings" -> dress, "AngleDefs" -> spAngles] ];
  Print[">>> ", tag, ": diagrams = ", Length[ntk[[1]]["Diagrams"]],
        "  NumTrace = ", Round[tA, 0.01], " s   MakeNTKernel(contract+compile+run) = ", Round[tG, 0.1], " s",
        "   TOTAL = ", Round[tA + tG, 0.1], " s"];
  {tA, tG}];

{aOn,  gOn}  = timeGen[True,  "COLLECTION ON ", "za3on"];
{aOff, gOff} = timeGen[False, "COLLECTION OFF", "za3off"];
Print[">>> speedup: NumTrace ", Round[aOff/aOn, 0.1], "x   MakeNTKernel ", Round[gOff/gOn, 0.1],
      "x   TOTAL ", Round[(aOff + gOff)/(aOn + gOn), 0.1], "x"];
