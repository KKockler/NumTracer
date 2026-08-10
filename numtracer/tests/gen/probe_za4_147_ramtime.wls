#!/usr/bin/env wolframscript
(* Stage 4 (b)-target: collected vs distributed RAM/time/trace-count on the DENSE za4_147 quark box
   (4 internal AqbqDirect147 vertices). Caps diagrams via NDIAG (env) so we time a subset, not an hour. *)
If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$tmp = "/tmp/za4rt"; Quiet@CreateDirectory[$tmp];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd147_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];
$RecursionLimit = 1000000;
SetEnvironment["NT_GEN_PROFILE" -> "2"];  (* per-wave RSS + phase A/B split from the generator *)

frame = NumTracer`sp4Frame[Global`p, Global`l1, Global`cos1, Global`cos2, Global`phi,
          Global`p1, Global`p2, Global`p3, Global`p4, Global`l1];
args  = {Global`l1, Global`cos1, Global`cos2, Global`phi, Global`p, Global`k};
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4, ZAqbq7, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};
spAngles4 = {Global`cosl1p1 -> Sqrt[1 - Global`cos1^2] Global`cos2,
  Global`cosl1p2 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 - 2 Sqrt[2 - 2 Global`cos2^2] Cos[Global`phi]),
  Global`cosl1p3 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 + Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] - Sqrt[3] Sin[Global`phi])),
  Global`cosl1p4 -> -(1/3) Sqrt[1 - Global`cos1^2] (Global`cos2 + Sqrt[2 - 2 Global`cos2^2] (Cos[Global`phi] + Sqrt[3] Sin[Global`phi]))};

(* Full permutation group of the identical external legs; the ClassTrans/transAAqbq
   projectors are covariant under it. Symmetries must be stated by hand since FunKit
   0bacdd1 -- without this the derivative yields every leg ordering as its own diagram
   (measured on za4_147: 6183 diagrams / 164 s, vs 294 / 7.6 s with it). *)
symLegs = FMakeSymmetryList[{A[i1], A[i2], A[i3], A[i4]}];
fRG  = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3], A[i4]}, "Symmetries" -> symLegs] // FTruncate // FRoute;
proj = FTerm[TBGetProjector["AAAAClassTrans", 1, {i1, i2, i3, i4} /. fRG["1-Loop"]["ExternalIndices"]]];
traceA4 = proj ** (fRG["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
flow = SPParam[traceA4 // dressingRules];

nKeep = With[{e = Environment["NDIAG"]}, If[StringQ[e] && IntegerQ[Quiet@ToExpression[e]], ToExpression[e], 4]];
(* cut the SAME PHYSICS pre-FromFunKit: keep nKeep loop terms; both modes then compute these, collected
   eagerly (few nets) vs distributed (each vertex sum expanded → many nets). *)
flowCut = If[nKeep > 0 && Head[flow] === Plus && Length[flow] > nKeep,
   Plus @@ Take[Sort[List @@ flow], nKeep], flow];
Print[">>> flow head=", Head[flow], " terms=", If[Head[flow] === Plus, Length[flow], 1], " -> cut ", If[Head[flowCut] === Plus, Length[flowCut], 1]];

run[coll_, ns_] := Module[{net, tA, ntk, tG},
  net = FromFunKit[flowCut, "DressingCollection" -> coll];
  {tA, ntk} = AbsoluteTiming[NumTrace[net, "Frame" -> frame, "Args" -> args, "DressingCollection" -> coll]];
  tG = First @ AbsoluteTiming[ MakeNTKernel[ntk,
     FileNameJoin[{$tmp, "gen_" <> ns <> ".cpp"}], FileNameJoin[{$tmp, ns <> "_k.hh"}],
     FileNameJoin[{$tmp, ns <> "_ks.hh"}], "Name" -> ns <> "_k", "Namespace" -> ns,
     "Dressings" -> dress, "AngleDefs" -> spAngles4] ];
  Print[">>> ", If[coll, "COLLECTED ", "DISTRIBUTED"], "  net terms=", If[Head[net] === Plus, Length[net], 1],
        "  diagrams=", Length[ntk[[1]]["Diagrams"]], "  NumTrace=", Round[tA, 0.01], "s  MakeNTKernel=", Round[tG, 0.1], "s"];
  {tA, tG}];

Print[">>> NDIAG cap = ", nKeep];
{aC, gC} = run[True, "za4c"];
{aD, gD} = run[False, "za4d"];
Print[">>> total: COLLECTED ", Round[aC + gC, 0.1], "s   DISTRIBUTED ", Round[aD + gD, 0.1], "s"];
