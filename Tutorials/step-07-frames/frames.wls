#!/usr/bin/env wolframscript
(* step-07 — Frames and kinematics.
 *
 * A frame is the choice of concrete components for every momentum in a network. This script
 *
 *   (a) PRINTS the component tables three of the builders produce, so you can see that a frame is
 *       nothing more mysterious than an Association  momentum-tag -> {c0, c1, c2, c3};
 *   (b) generates the SAME network twice — once on the symmetric-point 3-point frame `sp3Frame`,
 *       once on the general 3-point frame `gen3Frame` — so that the accompanying C++ checker can
 *       verify the documented relationship between them:
 *
 *           gen3Frame at  p1m == p2m,  cosP == -1/2   ==   sp3Frame.
 *
 * That identity is worth a test rather than a footnote: `gen3Frame` exists precisely because some
 * tensor bases are DEGENERATE at the symmetric point (their Gram determinant vanishes there), and
 * knowing the general frame reduces correctly is what lets you trust it off that slice.
 *
 * Needs a Wolfram kernel and FunKit. No TensorBases, no FORM, no DiFfRG.
 *)

If[Quiet[Needs["NumTracer`"]] === $Failed,
  Module[{dir = Environment["NUMTRACER_DIR"]},
    If[dir === $Failed, Print["step-07: set NUMTRACER_DIR"]; Exit[1]];
    Get[FileNameJoin[{dir, "mathematica", "NumTracer.m"}]]]];
Global`Nc = 3;

$out = Environment["NT_TUTORIAL_OUT"];
If[$out === $Failed, $out = Directory[]];

{p, l1, cos1, cos2} = {Global`p, Global`l1, Global`cos1, Global`cos2};
{p1m, p2m, cosP}    = {Global`p1m, Global`p2m, Global`cosP};
{q1, q2, q3, ql}    = {Global`q1, Global`q2, Global`q3, Global`ql};

(* ---------------------------------------------------------------------------------------------
 * @snip begin: show
 * A frame is DATA. Printing one is the fastest way to understand what a builder did — and the
 * fastest way to catch a frame that is not the kinematics you meant.
 *)
show[label_, fr_] := (
  Print["--- ", label, " ---"];
  Do[Print["  ", k, " = ", fr[k]], {k, Keys[fr]}]);

show["propFrame[p, l1, cos1, q1, ql]  (1 external, 1 angle)",
     propFrame[p, l1, cos1, q1, ql]];

show["sp3Frame[p, l1, cos1, cos2, q1, q2, q3, ql]  (3 externals at the symmetric point, 2 angles)",
     sp3Frame[p, l1, cos1, cos2, q1, q2, q3, ql]];

show["gen3Frame[p1m, p2m, cosP, l1, cos1, cos2, q1, q2, q3, ql]  (3 externals, general)",
     gen3Frame[p1m, p2m, cosP, l1, cos1, cos2, q1, q2, q3, ql]];
(* @snip end: show *)

(* ---------------------------------------------------------------------------------------------
 * @snip begin: nets
 * One network, two frames. The network mentions only momentum TAGS (q1, q2, ql) and index labels,
 * so it is entirely frame-independent — which is the point. Swapping the frame changes the
 * kinematics the kernel is evaluated at, not the algebra.
 *
 * The network is a 3-point-ish contraction: two different external legs tied through the loop
 * momentum's transverse projector, plus a scalar product so that the answer depends on the
 * external opening angle as well as on the loop angles.
 *
 * NOTE THE PARENTHESES. In a .wls script a line that is already a complete expression ENDS there,
 * so writing
 *     net = termA
 *         + termB;
 * silently assigns only termA and leaves `+termB` as a separate, discarded expression. Wrapping the
 * whole sum in ( ) is the fix. The "expected 2 diagrams" assertion below is what catches it.
 *)
net = ( ntVec[q1, mu] ntTransProj[ql, mu, nu] ntVec[q2, nu]
      + ntSP[q1, q2] ntVec[q3, rho] ntTransProj[ql, rho, sig] ntVec[q3, sig] );

frameSP  = sp3Frame[p, l1, cos1, cos2, q1, q2, q3, ql];
argsSP   = {l1, cos1, cos2, p};

frameGEN = gen3Frame[p1m, p2m, cosP, l1, cos1, cos2, q1, q2, q3, ql];
argsGEN  = {l1, cos1, cos2, p1m, p2m, cosP};
(* @snip end: nets *)

emit[fr_, ar_, nm_, ns_] := Module[{ntk, nd},
  ntk = NumTrace[net, "Frame" -> fr, "Args" -> ar];
  nd  = Length[ntk[[1]]["Diagrams"]];
  Print["step-07: ", nm, " diagrams = ", nd];
  (* Assert the EXPECTED count, not merely a nonzero one: the network above is a two-term sum, and
     a silently truncated sum (see the note on parentheses) would still give a perfectly valid
     kernel for the wrong network. *)
  If[nd =!= 2, Print["step-07 FAIL: expected 2 diagrams, got ", nd]; Exit[1]];
  MakeNTKernel[ntk,
    FileNameJoin[{$out, "gen_" <> ns <> ".cpp"}],
    FileNameJoin[{$out, nm <> ".hh"}],
    FileNameJoin[{$out, nm <> "s.hh"}],
    "Name" -> nm, "Namespace" -> ns, "Dressings" -> {}]];

emit[frameSP,  argsSP,  "sp3_kernel", "step07sp"];
emit[frameGEN, argsGEN, "gen3_kernel", "step07gen"];

Print["step-07: wrote sp3_kernel.hh / gen3_kernel.hh to ", $out];
