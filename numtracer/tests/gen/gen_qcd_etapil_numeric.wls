#!/usr/bin/env wolframscript
(* QCD pion anomalous dimension etaPiL (pion two-point wave-function flow) via the NUMERIC
   matrix-product backend. The pion self-energy gets its momentum (p^2) dependence purely from the
   quark loop (two pion-quark Yukawa vertices + the regularised quark propagator with a constituent
   mass Mq ~ hSigL*sqrt(rhoL/Nf)); the pointlike-meson loops are momentum-independent and drop out
   of the wave-function extraction (f(p)-f(0))/p^2 done before tracing. Emits the Numeric kernel
   (validated against the FORM oracle from gen_qcd_etapil_form.wls). Port of QCD.nb cells 14-16. Run from numtracer/:
     wolfram -script tests/gen/gen_qcd_etapil_numeric.wls                                          *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];

$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_meson_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetOptions[MakeNTKernel, Sequence @@ $ntShimOpts];

scalarParams = {etaPiL, etaSigL, d1V, d2V, rhoL};

(* frames (QCD.nb cell 7): NumTracer`propFrame head, p baked to 10^-3 for the wave-function flow *)
frame1 = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
frame0 = frame1 /. Global`p -> 10^(-3);

(* ---- pion two-point loop + wave-function extraction (QCD.nb cells 14-15) ----------------------- *)
fRGPiPi = FTakeDerivatives[WetterichEquation, {ΠL[i1], ΠL[i2]}] // FTruncate // FRoute;
traceΠL = FTerm[deltaAdjFlav[f1, f2], fRGPiPi["1-Loop"]["Expression"] /. diagRules] //. p2 -> -p1;
traceΠL = (Times @@ #1 &) /@ List @@ traceΠL;
traceΠL = Simplify[UseLorentzLinearity[(traceΠL - (traceΠL //. p1 -> 0))/sp[p1, p1]]];
traceΠL = FEx @@ FTerm /@ traceΠL;
Print["traceΠL terms = ", If[Head[traceΠL] === FEx, Length[traceΠL], 1]];

flow = PropParam[dressingRules[FTerm[-1] ** traceΠL]] /. p -> 10^(-3);
net  = FromFunKit[flow];
Print["etaPiL net top terms = ", If[Head[net] === Plus, Length[net], 1]];

ntk = NumTrace[net, "Frame" -> frame0, "Args" -> args0];

(* NUMERIC backend: signature kernel(l1, cos1, k, etaPiL, etaSigL, d1V, d2V, rhoL, <14 dressings>) *)
tN = First @ AbsoluteTiming[ MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_qcd_etapil_num.cpp"}],
  FileNameJoin[{$gen, "EtaPiL_num_kernel.hh"}],
  FileNameJoin[{$gen, "EtaPiL_num_kernels.hh"}],
  "Name" -> "EtaPiL_num_kernel", "Namespace" -> "etapil_num",
  "Dressings" -> $ntDressings, "ScalarParams" -> scalarParams] ];
Print["[time] Numeric backend (contract+compile+run): ", tN, " s"];
Print["etaPiL numeric kernels generated."];
