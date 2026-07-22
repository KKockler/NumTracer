#!/usr/bin/env wolframscript
(* Stage 4 north-star probe: FULL-BASIS ZAAqbq (AqbqDirect147 internal quark-gluon vertex).
   Generates collection ON vs OFF on the cut ZAAqbq1(small) topology and reports diagram count,
   generator source size, compile/run time, and phase-A table/RSS. Writes to /tmp so no committed
   fixture is touched. Run from numtracer/:  wolfram -script tests/gen/probe_zaaqbq147.wls *)
If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$tmp = "/tmp/zaaqbq147"; Quiet@CreateDirectory[$tmp];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_aaqbq147_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];
$RecursionLimit = 1000000;

frame = NumTracer`sp4Frame[Global`p, Global`l1, Global`cos1, Global`cos2, Global`phi,
          Global`p1, Global`p2, Global`p3, Global`p4, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`phi, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4, ZAqbq7, ZAAqbq1, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};
spAngles4 = {
  Global`cosl1p1 -> Sqrt[1 - Global`cos1^2] Global`cos2,
  Global`cosl1p2 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 - 2 Sqrt[2 - 2 Global`cos2^2] Cos[Global`phi]),
  Global`cosl1p3 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 + Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] - Sqrt[3] Sin[Global`phi])),
  Global`cosl1p4 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 + Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] + Sqrt[3] Sin[Global`phi]))};

fRGAAqbq = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], qb[i3], q[i4]}] // FTruncate // FRoute;
proj = FTerm[TBGetProjector["transAAqbqR", 1, {i1, i2, i4, i3} /. fRGAAqbq["1-Loop"]["ExternalIndices"]]];
trace = FTerm[proj, fRGAAqbq["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]];
flow = SPParam[dressingRules[trace]];

nKeep = With[{e = Environment["NT_TEST_NDIAG"]}, If[StringQ[e] && IntegerQ[Quiet@ToExpression[e]], ToExpression[e], 6]];

srcBytes[ns_] := Total[FileByteCount /@ FileNames[ns <> "_num*.cpp", $tmp]];

run[coll_, tag_, ns_] := Module[{net, tA, ntk, tG},
  net = FromFunKit[flow, "DressingCollection" -> coll];
  If[nKeep > 0 && Head[net] === Plus && Length[net] > nKeep, net = Plus @@ Take[Sort[List @@ net], nKeep]];
  {tA, ntk} = AbsoluteTiming[NumTrace[net, "Frame" -> frame, "Args" -> args, "DressingCollection" -> coll]];
  tG = First @ AbsoluteTiming[ MakeNTKernel[ntk,
     FileNameJoin[{$tmp, ns <> "_num.cpp"}], FileNameJoin[{$tmp, ns <> "_num_kernel.hh"}],
     FileNameJoin[{$tmp, ns <> "_num_kernels.hh"}],
     "Name" -> ns <> "_num_kernel", "Namespace" -> ns, "Dressings" -> dress, "AngleDefs" -> spAngles4] ];
  Print[">>> ", tag, ": diagrams = ", Length[ntk[[1]]["Diagrams"]],
        "  analyse = ", Round[tA, 0.01], " s  gen(contract+compile+run) = ", Round[tG, 0.1], " s",
        "  source = ", Round[srcBytes[ns]/1048576., 0.01], " MB  kernel = ",
        Round[FileByteCount[FileNameJoin[{$tmp, ns <> "_num_kernels.hh"}]]/1024., 0.1], " KB"];
  {tA, tG}];

runs = With[{e = Environment["NT_PROBE_RUNS"]}, If[StringQ[e], e, "both"]]; (* "on" | "off" | "both" *)
Print["=== FULL-BASIS ZAAqbq (AqbqDirect147 internal vertex), cut to ", nKeep, " diagrams  [runs=", runs, "] ==="];
If[runs =!= "off", run[True,  "COLLECTION ON ", "zaaqbq147on"]];
If[runs =!= "on",  run[False, "COLLECTION OFF", "zaaqbq147off"]];
Print["=== done ==="];
