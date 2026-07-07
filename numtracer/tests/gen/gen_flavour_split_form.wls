#!/usr/bin/env wolframscript
(* FormTracer (FORM) oracle for the QCD flavour-split scalar self-energy flow (Flavour_split) — the
   independent cross-check for the NUMERIC matrix-product kernel, re-anchoring the value check of
   tests/refshim/compare_flavour_split.cpp off the (removed) Dense DTensor oracle.

   Traces the SAME two-flavour quark loop that gen_flavour_split.wls builds — a scalar (sigma)
   two-point self-energy summed over two EXPLICIT flavours u and d, each its own diagram carrying
   its OWN symbolic runtime dressing (Gu vs Gd) and regulator-dot insertion (GuDot vs GdDot), only
   colour SU(3) folded — but with FORM (FormTrace) instead of the numeric backend. The per-flavour
   dressings Gu/Gd/GuDot/GdDot are kept SYMBOLIC so MakeKernel emits them as the same four runtime
   dressing parameters, making the FORM kernel share the numeric kernel's signature exactly:
       kernel(l1, cos1, p, Gu, Gd, GuDot, GdDot)
   Emits Flavour_split_form_kernel.hh into refshim/. Run from numtracer/:
     wolfram -script tests/gen/gen_flavour_split_form.wls                                           *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];

(* clean cold run: empty this kernel's trace cache so FormTrace reruns tform; TBCache stays warm. *)
$cacheDir = FileNameJoin[{$repo, "numtracer", "tests", "gen", "TraceCache_flavour_split_form"}];
If[DirectoryQ[$cacheDir], DeleteDirectory[$cacheDir, DeleteContents -> True]];
FSetCacheDirectory[$cacheDir];

(* ---- field space + truncation (identical to gen_flavour_split.wls) ------------------------------
   one scalar sigma + two quark flavours (Dirac index di, fundamental-colour index ci). *)
fields = <|"Commuting" -> {sig[p]},
  "Grassmann" -> {{qub[p, {di, ci}], qu[p, {di, ci}]}, {qdb[p, {di, ci}], qd[p, {di, ci}]}}|>;
truncation = <|
  GammaN -> {{sig, sig}, {sig, qub, qu}, {sig, qdb, qd}, {qub, qu}, {qdb, qd}},
  Propagator -> {{sig, sig}, {qub, qu}, {qdb, qd}},
  R -> {{qub, qu}, {qdb, qd}}, Rdot -> {{qub, qu}, {qdb, qd}}, Field -> {{}}|>;
FSetGlobalSetup[<|"FieldSpace" -> fields, "Truncation" -> truncation,
  "FeynmanRules" -> <|GammaN -> {}, Propagator -> {}, R -> {}, Rdot -> {}|>|>];

fRG = FTakeDerivatives[WetterichEquation, {sig[i1], sig[i2]}] // FTruncate // FRoute;
expr = fRG["1-Loop"]["Expression"];

(* hand-rolled diagrammatic rules, one set per flavour (identical to gen_flavour_split.wls). The
   ONLY per-flavour difference is the dressing name on the propagator / Rdot insertion. Yukawa:
   i gamma5 (Dirac) x colour delta. Quark line: massless q-slash x colour delta x per-flavour
   dressing.  The scalar sigma leg carries a momentum but NO index list, hence the leg pattern {_}. *)
yuk[fb_, f_] := GammaN[{sig, fb, f}, {{_}, {_, {dqb_, cqb_}}, {_, {dq_, cq_}}}] :>
  I gamma5[dqb, dq] deltaFundCol[cqb, cq];
prop[fb_, f_] := Propagator[{f, fb} | {fb, f}, {{m1_, {d1_, c1_}}, {m2_, {d2_, c2_}}}] :>
  With[{mu = Unique["lm"]},
    deltaFundCol[c1, c2] (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[GammaN, {fb, f}, 1, {m1, m2}]];
rdot[fb_, f_] := Rdot[{fb, f}, {{m1_, {d1_, c1_}}, {m2_, {d2_, c2_}}}] :>
  With[{mu = Unique["lm"]},
    deltaFundCol[c1, c2] (-I) gamma[mu, d1, d2] vec[m2, mu] * dressing[Rdot, {fb, f}, 1, {m1, m2}]];
myRules = {yuk[qub, qu], yuk[qdb, qd], prop[qub, qu], prop[qdb, qd], rdot[qub, qu], rdot[qdb, qd]};

(* per-flavour SYMBOLIC dressings -> the runtime kernel params Gu/Gd/GuDot/GdDot. *)
dressed = (expr /. myRules) //. {
  dressing[GammaN, {qub, qu}, 1, {m1_, m2_}] :> Gu[Sqrt[sp[m2, m2]]],
  dressing[GammaN, {qdb, qd}, 1, {m1_, m2_}] :> Gd[Sqrt[sp[m2, m2]]],
  dressing[Rdot,   {qub, qu}, 1, {m1_, m2_}] :> GuDot[Sqrt[sp[m2, m2]]],
  dressing[Rdot,   {qdb, qd}, 1, {m1_, m2_}] :> GdDot[Sqrt[sp[m2, m2]]]};
dressed = UseLorentzLinearity[dressed] //. {lf1 -> l1};

(* colour SU(3) baked; NO flavour group (no flavour heads). Nf enters only via the closed colour
   trace tr(1)=Nc per loop, which FORM folds numerically. *)
SetNc[3];

(* vacuum 1-angle scalar reduction: sp/cos -> {l1, p, cos1}, matching propFrame in the numeric
   script (copied from qcd_setup.m's PropParam). *)
PropParam[e_] := UseLorentzLinearity[e] //. {
  lf1 -> l1,
  sp[p1, p1] -> p^2, sp[l1, l1] -> l1^2,
  sp[l1, p1] -> l1 p cos[l1, p], sp[p1, l1] -> l1 p cos[l1, p],
  Sqrt[a_^2] :> a, (a_^2)^(n_/2) :> a^n,
  cos[l1, p] :> cos1};

$Assumptions = p > 0 && l1 > 0 && -1 < cos1 < 1;

{tTrace, Flow} = AbsoluteTiming[
  (FormTrace["Flavour_split_form", dressed] // PropParam // Simplify)];
Print["[time] FORM trace (FTakeDerivatives..FormTrace..PropParam): ", tTrace, " s"];

(* DiFfRG kernel parameter list matching the NUMERIC Flavour_split signature
   (Flavour_split_num_kernel.hh): kernel(l1, cos1, p, Gu, Gd, GuDot, GdDot). l1, cos1 are the
   integration variables; p is a plain double; Gu/Gd/GuDot/GdDot are shim SplineInterpolator1D
   references, in that order. *)
interpolatorType = "SplineInterpolator1D<double, LogarithmicCoordinates1D<double>, GPU_memory>";
splitDressings = {Gu, Gd, GuDot, GdDot};
splitParameterList = Join[
  {<|"Name" -> "p", "Type" -> "double"|>},
  Table[<|"Name" -> ToString[d], "Type" -> interpolatorType, "Const" -> True, "Reference" -> True|>, {d, splitDressings}]];

tKernel = First @ AbsoluteTiming[
  MakeKernel[Flow,
    "Name" -> "Flavour_split_form", "Integrator" -> "Integrator_p2_1ang", "d" -> 4, "AD" -> False,
    "ctype" -> "double", "Device" -> "GPU", "Type" -> "double",
    "Parameters" -> splitParameterList,
    "IntegrationVariables" -> {"l1", "cos1"},
    "Coordinates" -> {}, "CoordinateArguments" -> {}]];
Print["[time] FormTracer MakeKernel (FORM codegen + CSE): ", tKernel, " s"];

(* place the generated kernel header where the compare test includes it: swap DiFfRG physics
   includes for the standalone shim, matching the other refshim/*_form_kernel.hh oracles. *)
$src = FileNameJoin[{$repo, "flows", "Flavour_split_form", "kernel.hh"}];
$out = FileNameJoin[{$repo, "numtracer", "tests", "refshim", "Flavour_split_form_kernel.hh"}];
If[FileExistsQ[$src],
  Module[{txt},
    txt = Import[$src, "Text"];
    txt = StringReplace[txt,
      "#include \"DiFfRG/physics/interpolation.hh\"\n#include \"DiFfRG/physics/physics.hh\"" ->
        "#include \"shim.hpp\""];
    Export[$out, txt, "Text"];
    Print["wrote oracle: ", $out]],
  Print["WARN: ", $src, " not found"]];
Print["DONE Flavour_split FORM kernel."];
