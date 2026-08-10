#!/usr/bin/env wolframscript
(* step-11 — Projectors and tensor bases: assert, don't hope.
 *
 * step-10 projected the quark two-point flow onto element 1 of the "qbq" basis, using a projector
 * with a reversed index order and a leading minus sign — three conventions, none of them checkable
 * by NumTracer. This script shows the check that DOES catch a wrong projector, and it needs no
 * tracer, no flow, and no oracle: the DEFINING property of a dual basis,
 *
 *     <P_i, V_j>  =  delta_ij      for every pair (i, j)
 *
 * where P_i = TBGetProjector[...] and V_j = TBGetVertex[...]. If that fails, everything downstream
 * is wrong, and it fails LOUDLY at the top of your script instead of quietly in the physics.
 *
 * TIER C. Needs FunKit, TensorBases, DiFfRG and FORM. Not wired into ctest; run it as
 *
 *     NT_TUTORIAL_OUT=/tmp/out wolframscript -script orthonormality.wl
 *)

If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];

$repo = AbsoluteFileName@FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
SetDirectory[$repo];
DefineFormExecutable["/usr/bin/tform -w16"];
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];

$out = Environment["NT_TUTORIAL_OUT"];
If[$out === $Failed, $out = Directory[]];

(* --------------------------------------------------------------------------------------------
 * @snip begin: guard
 * FIRST, guard against a VACUOUS pass. If a basis never resolved — no cached Gram matrix and FORM
 * did not build one — TBGetProjector stays UNEVALUATED. It is then a symbolic head that happily
 * flows through the whole pipeline, and every check below "passes" while testing nothing.
 *
 * This is the single most important line in the script. A green orthonormality gate on an
 * unresolved basis is worse than no gate.
 *)
(* The two legs of a two-point function. Momentum conservation forces the second to be -p1, and the
   MOMENTUM TRAVELS WITH ITS LEG: each entry is one bundle {momentum, Dirac, colour, flavour}.
   idxP is idxV with the two bundles SWAPPED — see the demonstration further down for why. *)
idxV = {{Global`p1, Global`d1, Global`A1, Global`F1},
        {-Global`p1, Global`d2, Global`A2, Global`F2}};
idxP = Reverse[idxV];

If[! FreeQ[TBGetProjector["qbq", 1, idxP], _TBGetProjector] ||
   ! FreeQ[TBGetVertex["qbq", 1, idxV], _TBGetVertex],
  Print["step-11 FAIL: the \"qbq\" basis did not resolve — TBGetProjector/TBGetVertex are still ",
        "unevaluated heads. Every check below would pass vacuously. Fix the basis setup (TBCache ",
        "location, FORM executable) before trusting anything downstream."];
  Exit[1]];
Print["step-11: basis \"qbq\" resolved."];
(* @snip end: guard *)

(* --------------------------------------------------------------------------------------------
 * @snip begin: pairs
 * The check itself. "qbq" has two elements — the wave-function structure I*pdash and the mass
 * structure — so there are 2x2 = 4 pairs. FormTrace closes the Dirac/colour/flavour indices; the
 * result must be the identity matrix.
 *
 * THREE CONVENTIONS ARE BEING PINNED HERE AT ONCE, and each has cost somebody a week:
 *
 *  (a) A projector is dual to the VERTEX, not to the basis element. Use TBGetVertex — do NOT route
 *      this through TBGetInnerProduct, which for a vertex basis selects the canonical product used
 *      to project onto vertices and so mixes two conventions.
 *  (b) The index order of the projector is REVERSED relative to the vertex (idxP vs idxV above).
 *  (c) THE MOMENTUM TRAVELS WITH ITS LEG. The whole {momentum, Dirac, colour, flavour} bundle
 *      moves together. Swapping only the (d,A,F) triples while pinning the momenta transposes the
 *      Dirac structure instead — which shows up as sign flips on the C-odd structures.
 *)
n = 2;
gram[iP_, iV_] := Table[
   Simplify[FormTrace[TBGetProjector["qbq", i, iP] * TBGetVertex["qbq", j, iV]]],
   {i, n}, {j, n}];

mat = gram[idxP, idxV];
Print["step-11: <P_i, V_j> = ", InputForm[mat]];

If[mat =!= IdentityMatrix[n],
  Print["step-11 FAIL: <P_i, V_j> is not the identity. A wrong index order, a wrong dual, or a ",
        "momentum that did not travel with its leg. Everything downstream of this is wrong."];
  Exit[1]];
Print["step-11: orthonormality holds exactly."];
(* @snip end: pairs *)

(* --------------------------------------------------------------------------------------------
 * @snip begin: wrong
 * WHAT THE WRONG CONVENTION LOOKS LIKE. Give the projector the SAME leg order as the vertex — the
 * natural-looking thing to write — and the (1,1) entry comes back as
 *
 *     -gamma[p1, spi1, spi2]^2 / (4 sp[p1, p1])
 *
 * i.e. an UNCLOSED SPINOR LINE: the two gamma matrices never paired into a trace, so FormTrace left
 * a squared open-index object standing. That is the signature of a transposed Dirac structure, and
 * it is why step-10's projector carries {i2, i1} rather than {i1, i2}.
 *
 * Recognising this failure shape is worth more than memorising the rule.
 *)
wrong = gram[idxV, idxV];
Print["step-11: with the projector's legs NOT reversed, <P_i, V_j> = ", InputForm[wrong]];
If[wrong === IdentityMatrix[n],
  Print["step-11 FAIL: the deliberately wrong convention also gave the identity — this basis is ",
        "leg-symmetric, so the demonstration above is vacuous for it."]; Exit[1]];
(* @snip end: wrong *)

(* --------------------------------------------------------------------------------------------
 * @snip begin: weighted
 * SCALING UP: the weighted-trace trick.
 *
 * A basis with n structures has n^2 pairs. Checking each with its own generated kernel means n^2
 * committed headers — 144 for the full 12-element quark-gluon vertex basis. Instead contract the
 * WHOLE matrix against a fixed weight matrix in ONE net:
 *
 *     W  =  Sum_ij w[i,j] <P_i, V_j>     must equal     Sum_i w[i,i] .
 *
 * The weights are DISTINCT PRIMES, so no plausible error cancels: any single wrong entry shifts W
 * by its own weight, and the diagonal and off-diagonal weights are disjoint. One kernel covers all
 * n^2 pairs.
 *
 * (The per-pair matrix above is what you want when W is wrong and you need to know WHICH entry
 * moved — a diagnostic, not the gate.)
 *)
wOf[i_, j_] := Prime[n*(i - 1) + j];
expect = Total[Table[wOf[i, i], {i, n}]];
weighted = Total[Flatten[Table[
   wOf[i, j]*Simplify[FormTrace[TBGetProjector["qbq", i, idxP]*TBGetVertex["qbq", j, idxV]]],
   {i, n}, {j, n}]]];

Print["step-11: weighted trace = ", weighted, "   (expected ", expect, ")"];
If[weighted =!= expect,
  Print["step-11 FAIL: weighted trace mismatch."]; Exit[1]];
(* @snip end: weighted *)

Print["step-11: all projector checks passed."];
