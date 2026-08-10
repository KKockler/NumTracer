#!/usr/bin/env wolframscript
(* P3 — run the label census against the REAL pure-YM ZA4 net (the measured-broken flow).
   Emits no kernel; it only builds the net and runs NumTrace, which is where checkLabels fires. *)
If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "ym_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Print["FunKit = ", FunKit`$FunKitDirectory];

frame = NumTracer`sp4Frame[Global`p, Global`l1, Global`cos1, Global`cos2, Global`phi,
          Global`p1, Global`p2, Global`p3, Global`p4, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`phi, Global`p, Global`k};

(* Full permutation group of the identical external legs; the ClassTrans/transAAqbq
   projectors are covariant under it. Symmetries must be stated by hand since FunKit
   0bacdd1 -- without this the derivative yields every leg ordering as its own diagram
   (measured on za4_147: 6183 diagrams / 164 s, vs 294 / 7.6 s with it). *)
symLegs = FMakeSymmetryList[{A[i1], A[i2], A[i3], A[i4]}];
fRG  = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3], A[i4]}, "Symmetries" -> symLegs] // FTruncate // FPlot // FRoute // FPrint;
proj = FTerm[TBGetProjector["AAAAClassTrans", 1, {i1, i2, i3, i4} /. fRG["1-Loop"]["ExternalIndices"]]];
traceA4 = proj ** (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
net = FromFunKit[ SPParam[traceA4 // dressingRules] ];
Print["YM ZA4 net top terms = ", If[Head[net] === Plus, Length[net], 1]];

(* Raw evidence, independent of checkLabels: the global label tally over the WHOLE net.
   Expect noise here from un-distributed Plus vertices — this is only a cross-check. *)
tal = Sort[Tally[Flatten @ Cases[net, f_?NumTracer`Private`tensorQ :>
        NumTracer`Private`labelsOf[f], Infinity]], #1[[2]] > #2[[2]] &];
Print["top whole-net label counts (noisy by construction): ", Take[tal, UpTo[12]]];

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args];   (* checkLabels fires here *)
Print["CENSUS CLEAN: ", Length[ntk[[1]]["Diagrams"]], " diagrams, every label 1 or 2."];
