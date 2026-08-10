#!/usr/bin/env wolframscript
(* step-19 — Fixed-component gammas: gamma^0 is a slash against a unit basis vector.
 *
 * Four-quark Fierz bases are written in a finite-T 3+1 split, so their structures carry gamma^0 —
 * a Lorentz index pinned to a CONCRETE COMPONENT rather than a contraction label. The DSL had no
 * notion of that: `ntGamma[0, d1, d2]` was read as "gamma^mu with a label named 0", the literal 0
 * then occurred in many factors at once, and NumTrace aborted:
 *
 *     NumTrace::privclash: label(s) 0 are private dummies of one factor but also occur outside it
 *
 * The guard was behaving correctly — the FEATURE was missing. It is supplied by a rewrite inside
 * NumTrace (expandFixedComponents): a fixed component is a contraction with the constant unit basis
 * vector e_i, so gamma^i becomes an ordinary SLASH against e_i. No new token kind, and ZERO C++
 * change. This script pins that rewrite against analytic oracles.
 *
 * THE FRAME IS propFrameFT, whose TEMPORAL components are nonzero. At a T=0 vacuum frame every q_0
 * vanishes and half the oracles below would collapse to 0 — the gate would pass vacuously no matter
 * what the tracer did.
 *
 * Needs a Wolfram kernel and FunKit. No TensorBases, no FORM, no DiFfRG.
 *)

If[Quiet[Needs["NumTracer`"]] === $Failed,
  Module[{dir = Environment["NUMTRACER_DIR"]},
    If[dir === $Failed, Print["step-19: set NUMTRACER_DIR"]; Exit[1]];
    Get[FileNameJoin[{dir, "mathematica", "NumTracer.m"}]]]];
Global`Nc = 3;

$out = Environment["NT_TUTORIAL_OUT"];
If[$out === $Failed, $out = Directory[]];

{p0, p, l0, l1, cos1} = {Global`p0, Global`p, Global`l0, Global`l1, Global`cos1};
{q1, ql, e0}          = {Global`q1, Global`ql, Global`e0};
args = {p0, p, l0, l1, cos1};

(* q1 = {p0, p, 0, 0};   ql = {l0, l1 cos1, l1 sin1, 0} *)
frame = propFrameFT[p0, p, l0, l1, cos1, q1, ql];

(* @snip begin: refframe
 * THE REFERENCE FRAME. The same traces can be written WITHOUT the fixed-component feature at all,
 * by adding an explicit unit vector as an ordinary user momentum: gamma^0 == gamma^mu contracted
 * with (e_0)_mu. Grading the two spellings against EACH OTHER is convention-independent — it cannot
 * be fooled by a sign or metric-convention error that I also made in the hand-written oracle, which
 * is the failure mode that produces a false green.
 *)
frameRef = Append[frame, e0 -> {1, 0, 0, 0}];
(* @snip end: refframe *)

emit[net_, nm_, ns_, fr_] := Module[{ntk, nd},
  ntk = NumTrace[net, "Frame" -> fr, "Args" -> args, "Dressings" -> {}];
  nd  = Length[ntk[[1]]["Diagrams"]];
  If[nd =!= 1, Print["step-19 FAIL: ", nm, " gave ", nd, " diagrams, expected 1"]; Exit[1]];
  MakeNTKernel[ntk,
    FileNameJoin[{$out, "gen_" <> ns <> ".cpp"}],
    FileNameJoin[{$out, nm <> ".hh"}],
    FileNameJoin[{$out, nm <> "s.hh"}],
    "Name" -> nm, "Namespace" -> ns, "Dressings" -> {}];
  Print["step-19: emitted ", nm]];

(* @snip begin: rewrite
 * ASSERTION A — the rewrite happens at all, and puts a unit vector in the env.
 * Pre-fix this NumTrace call aborted outright. Post-fix the fixed index must be gone from the label
 * machinery AND the unit vector must have reached the env as a momentum: without an env entry,
 * the slash emission has nothing to reference.
 *)
netA = ntGamma[0, d1, d2] ntGamma[a1, d2, d1] ntVec[q1, a1];
ntkA = NumTrace[netA, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
If[! MemberQ[Keys[ntkA[[1]]["Env"]], ntUnitVec[0]],
  Print["step-19 FAIL: e_0 did not reach the env; a fixed-component gamma cannot be emitted as a ",
        "slash without one. Env = ", Keys[ntkA[[1]]["Env"]]]; Exit[1]];
Print["step-19: fixed component rewritten, ntUnitVec[0] present in the env."];
(* @snip end: rewrite *)

(* @snip begin: components
 * ASSERTION B — 0-based component indexing.
 *   tr(gamma^0 q1slash) = 4 q1_0 = 4 p0     and     tr(gamma^1 q1slash) = 4 q1_1 = 4 p .
 * In propFrameFT, p0 and p are INDEPENDENT kernel arguments, so an off-by-one in the unit vector
 * (e_i placed at slot i+1, or 1-based components) SWAPS the two answers and is caught. That is the
 * whole reason both are emitted rather than just one.
 *)
emit[ntGamma[0, d1, d2] ntGamma[a1, d2, d1] ntVec[q1, a1], "Fx_t_kernel",  "fx_t",  frame];
emit[ntGamma[1, d1, d2] ntGamma[a1, d2, d1] ntVec[q1, a1], "Fx_x_kernel",  "fx_x",  frame];

(* ASSERTION C — a fixed component contracted with ITSELF: tr(gamma^0 gamma^0) = 4 (e_0.e_0) = 4.
   Two fixed components on one chain get two DISTINCT fresh dummies but the SAME env momentum, so
   this pins that the rewrite neither fuses the dummies (making the pair a self-contraction) nor
   allocates e_0 twice. *)
emit[ntGamma[0, d1, d2] ntGamma[0, d2, d1], "Fx_tt_kernel", "fx_tt", frame];
(* @snip end: components *)

(* @snip begin: crosscheck
 * ASSERTION D/E — the convention-independent cross-check.
 * D: the full four-gamma trace with fixed components INTERLEAVED with ordinary slashes,
 *      tr(gamma^0 q1slash gamma^0 qlslash) = 4 [2 q1_0 ql_0 - q1.ql] = 4 [p0 l0 - p l1 cos1].
 *    Kinematics-dependent in every argument, so it cannot fold to a constant.
 * E: THE SAME TRACE, with gamma^0 written the long way as gamma^mu against an explicit unit-vector
 *    momentum supplied by the user's frame. This uses none of the new code path, so agreement
 *    between D and E proves the rewrite is faithful INDEPENDENTLY of whether my Euclidean sign
 *    conventions in the oracle above are right.
 *)
netD = ntGamma[0, d1, d2] ntGamma[a1, d2, d3] ntVec[q1, a1] *
       ntGamma[0, d3, d4] ntGamma[a2, d4, d1] ntVec[ql, a2];
emit[netD, "Fx_4_kernel", "fx_4", frame];

netE = ntGamma[b1, d1, d2] ntVec[e0, b1] ntGamma[a1, d2, d3] ntVec[q1, a1] *
       ntGamma[b2, d3, d4] ntVec[e0, b2] ntGamma[a2, d4, d1] ntVec[ql, a2];
emit[netE, "Fx_4ref_kernel", "fx_4ref", frameRef];
(* @snip end: crosscheck *)

(* @snip begin: metric
 * ASSERTION F — a fixed component on a METRIC, not a gamma. The rewrite is not gamma-specific:
 *   g^{0 nu} q1_nu = q1_0 = p0,  so  tr(gamma^0 qlslash) g^{0 nu} q1_nu = 4 ql_0 p0 = 4 l0 p0.
 *)
netF = ntMetric[0, nu] ntVec[q1, nu] *
       ntGamma[0, d1, d2] ntGamma[a1, d2, d1] ntVec[ql, a1];
emit[netF, "Fx_g_kernel", "fx_g", frame];
(* @snip end: metric *)

Print["step-19: wrote 6 fixed-component kernels to ", $out];
