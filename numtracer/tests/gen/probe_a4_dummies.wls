#!/usr/bin/env wolframscript
(* P0 — DIAGNOSTIC ONLY, no tracer, no FORM, no loop.
   Does TensorBases uniquify AAAAClassTrans's internal Lorentz dummies v1..v4?
   AAAAClassTrans (bases/AAAAClassTrans.m:60-66) uses v1..v4 as INTERNAL dummies, while
   AAAAClass (bases/AAAAClass.m:60-63) uses v1..v4 as its EXTERNAL leg names. Same global
   symbols. If fixIndices/GetClosedIndices fails to rename the projector's dummies, the two
   objects collide the moment they meet in one diagram — which is exactly and only the
   rank-4-projector + AAAA-vertex case that is measured broken. *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
SetDirectory[$repo];
Get[FileNameJoin[{$repo, "numtracer", "tests", "ym_setup.m"}]];

Print["FunKit      = ", FunKit`$FunKitDirectory];
Print["TensorBases = ", If[ValueQ[TensorBases`$TensorBasesDirectory], TensorBases`$TensorBasesDirectory, "?"]];
Print[];

(* heads/wrappers that legitimately appear in both and are not indices *)
$boring = {Plus, Times, Power, List, Rational, Complex, Integer, Real, Symbol,
           Global`transProj, Global`TBdeltaLorentz, Global`TBF, Global`TBvec, Global`TBsp,
           Global`color, Global`Nc, Global`p1, Global`p2, Global`p3, Global`p4};

report[nm_, a_, b_, externals_] := Module[{sa, sb, shared},
  sa = DeleteDuplicates@Cases[a, _Symbol, Infinity];
  sb = DeleteDuplicates@Cases[b, _Symbol, Infinity];
  shared = Complement[Intersection[sa, sb], Join[$boring, externals]];
  Print["--- ", nm];
  Print["    A symbols = ", Sort[Complement[sa, $boring]]];
  Print["    B symbols = ", Sort[Complement[sb, $boring]]];
  Print["    SHARED non-external symbols (MUST be {}): ", shared];
  Print[If[shared === {}, "    => clean", "    => *** COLLISION ***"]];
  shared];

idx4 = {Global`i1, Global`i2, Global`i3, Global`i4};
idx3 = {Global`i1, Global`i2, Global`i3};

projA4 = TBGetProjector["AAAAClassTrans", 1, idx4];
vertA4 = TBGetBasisElement["AAAAClass", 1, idx4];
projA3 = TBGetProjector["AAAClass", 1, idx3];
vertA3 = TBGetBasisElement["AAAClass", 1, idx3];

Print["===== P0: cross-object dummy collisions ====="];
Print[];
bad4 = report["rank-4 projector (AAAAClassTrans) vs AAAA vertex (AAAAClass)  <-- THE BROKEN CASE",
              projA4, vertA4, idx4];
Print[];
bad3 = report["rank-3 projector (AAAClass) vs AAA vertex (AAAClass)  <-- known-good control",
              projA3, vertA3, idx3];
Print[];
badX = report["rank-3 projector (AAAClass) vs AAAA vertex (AAAAClass)  <-- known-good control (dtZA3)",
              projA3, vertA4, DeleteDuplicates@Join[idx3, idx4]];
Print[];

(* Independent view: two INDEPENDENT instantiations of the same object must share no dummy. *)
Print["===== P0b: two independent instantiations of AAAAClassTrans ====="];
p1a = TBGetProjector["AAAAClassTrans", 1, idx4];
p1b = TBGetProjector["AAAAClassTrans", 1, idx4];
Print["    identical expressions? ", p1a === p1b, "  (True => dummies are NOT freshened per call)"];
report["AAAAClassTrans instance #1 vs instance #2", p1a, p1b, idx4];
Print[];

Print["===== VERDICT ====="];
Print["broken-case collision : ", bad4];
Print["control  collisions   : ", bad3, "  ", badX];
If[bad4 =!= {} && bad3 === {},
  Print[">>> CONFIRMED: the rank-4 projector leaks un-uniquified dummies ", bad4,
        " that collide with the AAAA vertex. Controls are clean. This is the root cause."],
  Print[">>> Not confirmed by P0 alone; proceed to P1/P2."]];
