#!/usr/bin/env wolframscript
(* step-08 — The DSL, head by head.
 *
 * Every network NumTracer contracts is built from a small vocabulary of tensor HEADS. This script
 * exercises them one at a time: each head gets a network small enough that its value is obvious,
 * and the accompanying C++ checker compares the generated kernel against that closed form.
 *
 * All of them share the one-angle frame
 *
 *     q1 = (p, 0, 0, 0),        ql = (l1 c, l1 s, 0, 0),   c = cos1, s = sqrt(1 - c^2)
 *
 * so the invariants are  q1.q1 = p^2,  ql.ql = l1^2,  q1.ql = p l1 c.
 *
 * Needs a Wolfram kernel and FunKit. No TensorBases, no FORM, no DiFfRG.
 *)

If[Quiet[Needs["NumTracer`"]] === $Failed,
  Module[{dir = Environment["NUMTRACER_DIR"]},
    If[dir === $Failed, Print["step-08: set NUMTRACER_DIR"]; Exit[1]];
    Get[FileNameJoin[{dir, "mathematica", "NumTracer.m"}]]]];
Global`Nc = 3;

$out = Environment["NT_TUTORIAL_OUT"];
If[$out === $Failed, $out = Directory[]];

{p, l1, cos1} = {Global`p, Global`l1, Global`cos1};
{q1, ql}      = {Global`q1, Global`ql};

frame = propFrame[p, l1, cos1, q1, ql];
args  = {l1, cos1, p};

(* Emit one micro-kernel. The diagram-count assertion is the cheap guard every script should carry:
   a network that analyses to ZERO diagrams has been pruned by a rule you did not intend (an odd
   Dirac trace, a vanishing colour fold), and emitting it would produce a kernel that silently
   returns nothing. *)
emit[net_, nm_, ns_, dress_ : {}] := Module[{ntk},
  ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> dress];
  If[Length[ntk[[1]]["Diagrams"]] < 1,
    Print["step-08 FAIL: ", nm, " analysed to 0 diagrams"]; Exit[1]];
  Print["step-08: ", nm, " diagrams = ", Length[ntk[[1]]["Diagrams"]]];
  MakeNTKernel[ntk,
    FileNameJoin[{$out, "gen_" <> ns <> ".cpp"}],
    FileNameJoin[{$out, nm <> ".hh"}],
    FileNameJoin[{$out, nm <> "s.hh"}],
    "Name" -> nm, "Namespace" -> ns, "Dressings" -> dress]];

(* ============================================================================================
 * @snip begin: lorentz
 * THE LORENTZ SECTOR — metrics, vectors, projectors. Index labels repeat => summed.
 *)

(* ntMetric[mu, nu] is the Euclidean delta_{mu nu}. Sandwiched between two vectors it is the plain
   scalar product:  q1_mu delta_{mu nu} q1_nu = q1.q1 = p^2. *)
emit[ntMetric[mu, nu] ntVec[q1, mu] ntVec[q1, nu], "H_metric_kernel", "h_metric"];

(* ntTransProj[q, mu, nu] = delta_{mu nu} - q_mu q_nu / q^2, the transverse projector on q.
   q1.P^T(ql).q1 = p^2 - (p l1 c)^2 / l1^2 = p^2 (1 - c^2). *)
emit[ntVec[q1, mu] ntTransProj[ql, mu, nu] ntVec[q1, nu], "H_transproj_kernel", "h_transproj"];

(* ntLongProj[q, mu, nu] = q_mu q_nu / q^2, the complement.  q1.P^L(ql).q1 = p^2 c^2.
   The two must add to the metric case above — the checker verifies exactly that. *)
emit[ntVec[q1, mu] ntLongProj[ql, mu, nu] ntVec[q1, nu], "H_longproj_kernel", "h_longproj"];
(* @snip end: lorentz *)

(* ============================================================================================
 * @snip begin: dirac
 * THE DIRAC SECTOR — spinor axes are named like Lorentz labels, and a chain closes into a trace
 * when the spinor labels form a cycle.
 *)

(* ntDeltaDirac[din, dout] is the spinor identity. Two of them close a loop: tr(1) = 4.
   A constant kernel, but a useful one — a wrong tr(1) silently rescales every flow by a power of 4. *)
emit[ntDeltaDirac[d1, d2] ntDeltaDirac[d2, d1], "H_diracid_kernel", "h_diracid"];

(* ntGamma[mu, din, dout] is gamma^mu carrying a Lorentz axis and two spinor axes. The spinor
   labels d1 -> d2 -> d1 close the trace; the Lorentz axes a, b are contracted with two momenta:
   tr(gamma^a gamma^b) q1_a ql_b = 4 (q1.ql) = 4 p l1 c. *)
emit[ntGamma[a, d1, d2] ntGamma[b, d2, d1] ntVec[q1, a] ntVec[ql, b],
     "H_dirac2_kernel", "h_dirac2"];

(* Four gammas: tr(g^a g^b g^c g^d) q1_a ql_b q1_c ql_d
     = 4[ (q1.ql)^2 - (q1.q1)(ql.ql) + (q1.ql)^2 ]
     = 4[ 2 p^2 l1^2 c^2 - p^2 l1^2 ].
   This is the first net whose value is not guessable at a glance, which is exactly why it is here. *)
emit[ntGamma[a, d1, d2] ntGamma[b, d2, d3] ntGamma[c, d3, d4] ntGamma[e, d4, d1] *
       ntVec[q1, a] ntVec[ql, b] ntVec[q1, c] ntVec[ql, e],
     "H_dirac4_kernel", "h_dirac4"];
(* @snip end: dirac *)

(* ============================================================================================
 * @snip begin: sun
 * THE SU(N) SECTOR — every head carries its RANK as the leading argument, so a colour SU(3) and a
 * flavour SU(2) sector can coexist in one network without interfering.
 *)

(* ntSUNT[N, a, i, j] = (T^a)_{ij}.  tr(T^a T^a) = (N^2-1)/2 = 4 for SU(3) — the C_F of step-02,
   now generated rather than hand-folded. *)
emit[ntSUNT[3, a, ii, jj] ntSUNT[3, a, jj, ii] ntSP[q1, q1], "H_sunT_kernel", "h_sunt"];

(* ntSUNf[N, a, b, c] = f^{abc}.  f^{abc} f^{abc} = N(N^2-1) = 24. *)
emit[ntSUNf[3, a, b, c] ntSUNf[3, a, b, c] ntSP[q1, q1], "H_sunf_kernel", "h_sunf"];
(* @snip end: sun *)

(* ============================================================================================
 * @snip begin: scalars
 * SCALARS — these never contract. They multiply a diagram and are emitted as ordinary arithmetic.
 *
 * ntSP[a, b] is a scalar product resolved by the frame.
 *
 * A DRESSING is just an ordinary Mathematica head applied to its arguments — `myZ[p]` — that you
 * additionally name in "Dressings". That declaration is what turns it into a kernel parameter
 * (`const auto& myZ`) and emits the call site as `myZ(p)`. The engine never looks inside one: it is
 * an opaque callable the consumer supplies, which is how propagator dressings, regulators and
 * interpolators enter a kernel.
 *
 * (There is also a head `ntDress[head, args...]` whose usage message promises the same thing. It
 * does NOT currently work — it emits a literal `ntDress(head, args)` call that no runtime header
 * declares, so the kernel fails to compile. Use the plain-symbol form below; it is also what every
 * production flow uses.)
 *)
emit[Global`myZ[p] ntSP[q1, ql] ntMetric[mu, nu] ntVec[q1, mu] ntVec[q1, nu],
     "H_dress_kernel", "h_dress", {Global`myZ}];
(* @snip end: scalars *)

Print["step-08: wrote 9 head kernels to ", $out];
