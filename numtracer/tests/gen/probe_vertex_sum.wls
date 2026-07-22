#!/usr/bin/env wolframscript
(* Probe: dump the raw quark-gluon vertex sum (AqbqDirect147 = T1+T4+T7) as it appears in the net,
   so Stage-4 CP3 can be written against ground truth. Based on gen_qcd_za3_147_numeric.wls. *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd147_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
$RecursionLimit = 1000000;

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
          Global`p1, Global`p2, Global`p3, Global`l1];
dress = {ZA3, ZAcbc, ZA4, ZAqbq1, ZAqbq4, ZAqbq7, dtZc, Zc, dtZA, ZA, dtZq, Zq, Mq};

fRGA3 = FTakeDerivatives[WetterichEquation, {A[i1], A[i2], A[i3]}] // FTruncate // FRoute;
projA3 = FTerm[TBGetProjector["AAAClassTrans", 1, {i1, i2, i3} /. fRGA3["1-Loop"]["ExternalIndices"]]];
traceA3 = projA3 ** (fRGA3["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
flow = SPParam[traceA3 // dressingRules];
(* capture the head-rewritten expression BEFORE expandBridges distributes it, by stubbing
   expandBridges/contractFlavour to Identity (OwnValue substitution beats their DownValues). *)
net = Block[{NumTracer`Private`expandBridges = Identity, NumTracer`Private`contractFlavour = Identity},
  FromFunKit[flow, "DressingCollection" -> True]];
Print["pre-distribution expr head = ", Head[net], "  top terms = ", If[Head[net] === Plus, Length[net], 1]];

(* Use NumTracer's own private classifiers to find vertex sums. *)
tensorQ = NumTracer`Private`tensorQ;
freeIdx = NumTracer`Private`freeIdx;
labelsOf = NumTracer`Private`labelsOf;
spinorLabelsHead = NumTracer`Private`spinorLabelsHead;
dressedStructureSumQ = NumTracer`Private`dressedStructureSumQ;
diracNumeratorSumQ = NumTracer`Private`diracNumeratorSumQ;

allSpinorLabels = NumTracer`Private`allSpinorLabels;
fundamentalSUNQ = NumTracer`Private`fundamentalSUNQ;
adjointSUNQ = NumTracer`Private`adjointSUNQ;

(* colour labels anywhere under p (on any SU(N) head) *)
colourLabels[p_] := DeleteDuplicates @ Flatten @ Cases[p,
  h_ /; (fundamentalSUNQ[h] || adjointSUNQ[h]) :> labelsOf[h], {0, Infinity}];
(* free LORENTZ legs = free indices that are neither spinor nor colour (the open gluon axis lives here) *)
lorentzFreeLegs[p_] := Complement[freeIdx[p], allSpinorLabels[p], colourLabels[p]];

diracHeadQ[p_] := ! FreeQ[p, _NumTracer`ntGamma | _NumTracer`ntGamma5 | _NumTracer`ntSigma | _NumTracer`ntDeltaDirac];
(* the true quark-gluon vertex: a dressed Dirac structure sum with an OPEN LORENTZ leg *)
vertexCandQ[p_Plus] := TrueQ[dressedStructureSumQ[p]] && diracHeadQ[p] && Length[lorentzFreeLegs[p]] >= 1;
vertexCandQ[_] := False;

cands = DeleteDuplicates @ Cases[net, p_Plus /; vertexCandQ[p] :> p, {0, Infinity}];
Print["\n#### vertex-sum (open-Lorentz-leg) candidates found = ", Length[cands]];
Do[
  Print["\n================ vertex candidate ", i, " (nterms=", Length[cands[[i]]], ") ================"];
  Print["freeIdx        = ", freeIdx[cands[[i]]]];
  Print["lorentzFreeLegs= ", lorentzFreeLegs[cands[[i]]]];
  Print["per-term:"];
  Do[Print["  term ", j, ": ", FullForm[(List @@ cands[[i]])[[j]]]], {j, 1, Length[cands[[i]]]}];
, {i, 1, Length[cands]}];
