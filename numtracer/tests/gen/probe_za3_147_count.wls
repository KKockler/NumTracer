#!/usr/bin/env wolframscript
(* Quantify the Stage 4 collection win on za3_147: #diagrams collected vs distributed (analyse only). *)
If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd147_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
$RecursionLimit = 1000000;
frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2, Global`p1, Global`p2, Global`p3, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4, ZAqbq7, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};
(* Full permutation group of the identical external legs; the ClassTrans/transAAqbq
   projectors are covariant under it. Symmetries must be stated by hand since FunKit
   0bacdd1 -- without this the derivative yields every leg ordering as its own diagram
   (measured on za4_147: 6183 diagrams / 164 s, vs 294 / 7.6 s with it). *)
symLegs = FMakeSymmetryList[{A[i1], A[i2], A[i3]}];
fRGA3 = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3]}, "Symmetries" -> symLegs] // FTruncate // FRoute;
projA3 = FTerm[TBGetProjector["AAAClassTrans", 1, {i1, i2, i3} /. fRGA3["1-Loop"]["ExternalIndices"]]];
traceA3 = projA3 ** (fRGA3["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
flow = SPParam[traceA3 // dressingRules];

countDiags[coll_] := Module[{net, ntk},
  net = FromFunKit[flow, "DressingCollection" -> coll];
  ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "DressingCollection" -> coll];
  {If[Head[net] === Plus, Length[net], 1], Length[ntk[[1]]["Diagrams"]]}];

{netOn, diagOn}   = countDiags[True];
{netOff, diagOff} = countDiags[False];
Print["COLLECTION ON : net top terms = ", netOn,  "  diagrams = ", diagOn];
Print["COLLECTION OFF: net top terms = ", netOff, "  diagrams = ", diagOff];
Print["diagram reduction factor = ", N[diagOff/diagOn, 3], "x"];
