#!/usr/bin/env wolframscript
(* Codegen gate for a MULTI-TERM (Plus) Dirac projector of MIXED gamma parity: a hand-built DSL net
   (no FunKit / no TensorBases / no FORM / no dressings) reproducing the bug where the odd-trace
   filter in DSL.m judged parity on a DIAGRAM-GLOBAL gamma count and so discarded the whole diagram.

   The projector Plus below carries ALL-NUMERIC coefficients and lives in a single sector, so
   distributeQ is False (expandBridges leaves it EAGER for et::add) and collectibleDiracSumQ is
   False (that predicate is gated on dressedStructureSumQ, i.e. a NON-numeric coefficient). It is
   therefore neither distributed nor exempted — exactly the crack the bug fell through.

     proj = 1 + q3_slash                     branch gamma-counts 0 and 1  -> MIXED parity
     d1   = tr[ proj . q1_slash . ql_slash ]
          = tr[q1_slash ql_slash] + tr[q3_slash q1_slash ql_slash]
          = 4 (q1.ql)             + 0                                    (odd trace vanishes)

   The diagram-global gamma count is 0/1 + 1 + 1 = 3 -> ODD, so the pre-fix Select dropped the
   ENTIRE diagram, taking the non-vanishing EVEN branch with it: 0 diagrams, an identically-zero
   kernel, and no error raised. The correct verdict is per-branch: drop only when EVERY branch is
   odd. compare_multiterm_proj_num.cpp checks the emitted kernel against 4 (q1.ql) analytically,
   so a kernel that is merely non-empty but wrong still fails.

   The net depends on the loop angle through ql, so the surviving branch is genuinely exercised
   over the frame rather than folding to a constant. *)
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
Global`Nc = 3; (* SetNc lives in TensorBases (not loaded here); NumTrace reads Global`Nc *)

frame = NumTracer`sp3Frame[Global`p, Global`l1, Global`cos1, Global`cos2,
                                   Global`q1, Global`q2, Global`q3, Global`ql];
args  = {Global`l1, Global`cos1, Global`cos2, Global`p};

(* trailing '+' keeps the statement open across the line break (a leading '+' would end it) *)
proj = ntDeltaDirac[d1, d2] +
       ntGamma[c1, d1, d2] ntVec[Global`q3, c1];
net  = proj * ntGamma[a1, d2, d3] ntVec[Global`q1, a1] *
              ntGamma[a2, d3, d1] ntVec[Global`ql, a2];

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];

(* Assertion A: the mixed-parity diagram must SURVIVE the odd-trace filter. Pre-fix this is 0. *)
If[Length[ntk[[1]]["Diagrams"]] < 1,
  Print["FAIL: the mixed-parity multi-term Dirac Plus was dropped by the odd-gamma filter ",
        "(0 diagrams). The filter must judge parity per BRANCH, not on the diagram-global ",
        "gamma count."]; Exit[1]];
Print["[gen] ok: mixed-parity multi-term projector survives the filter (",
      Length[ntk[[1]]["Diagrams"]], " diagram(s))"];

(* Assertion C (the other direction — the fix must not become a blanket exemption): a Plus whose
   branches are ALL odd once distributed is genuinely zero and must still be pruned.
     projZ = 1 + q1_slash q2_slash            branch gamma-counts 0 and 2
     dZ    = tr[ projZ . ql_slash ]           -> branch totals 1 and 3, BOTH odd -> identically 0
   (Its diagram-global count is also 3, so this case was pruned correctly by accident before; it is
   here to pin that the branch-aware rule keeps pruning it.) *)
projZ = ntDeltaDirac[e1, e2] +
        ntGamma[b1, e1, f1] ntVec[Global`q1, b1] ntGamma[b2, f1, e2] ntVec[Global`q2, b2];
netZ  = projZ * ntGamma[b3, e2, e1] ntVec[Global`ql, b3];
ntkZ  = NumTrace[netZ, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
If[Length[ntkZ[[1]]["Diagrams"]] =!= 0,
  Print["FAIL: an all-odd-branch diagram (identically zero) was NOT pruned (",
        Length[ntkZ[[1]]["Diagrams"]], " diagram(s)) — the odd-trace pruning has been ",
        "over-broadened into a blanket exemption."]; Exit[1]];
Print["[gen] ok: all-odd-branch diagram still pruned (0 diagrams)"];

(* Assertion D: POWERS must not blow up the parity fold. A scalar power like sp[q2,q3]^2 (ubiquitous
   in Gram determinants / projector normalisations) reaches diracParities as Power[b, n]. Folding it
   as `Times @@ ConstantArray[b, n]` recurses forever, because Times re-collapses the copies back to
   b^n — it hits $RecursionLimit and feeds garbage to Union, silently corrupting every verdict in the
   diagram. Same net as above, scaled by a squared scalar and a squared tensor self-contraction. *)
projP = ntDeltaDirac[g1, g2] + ntGamma[h1, g1, g2] ntVec[Global`q3, h1];
netP  = Global`p^2 * (1 + Global`p^2)^3 *
        projP * ntGamma[h2, g2, g3] ntVec[Global`q1, h2] *
                ntGamma[h3, g3, g1] ntVec[Global`ql, h3];
ntkP = Quiet@Check[NumTrace[netP, "Frame" -> frame, "Args" -> args, "Dressings" -> {}], $Failed];
If[ntkP === $Failed || Length[ntkP[[1]]["Diagrams"]] < 1,
  Print["FAIL: a diagram carrying Power[...] factors did not survive — the parity fold most likely ",
        "recursed on Power (Times re-collapses ConstantArray copies) and corrupted the verdict."];
  Exit[1]];
Print["[gen] ok: Power factors fold correctly (", Length[ntkP[[1]]["Diagrams"]], " diagram(s))"];

(* Assertion E: a COLLAPSED spinor loop must keep its tr(1) = 4. orderDiracFacs drops δ "connector"
   factors, so a loop built only from ntDeltaDirac arrives token-free; the runtime's split_loops
   discards empty segments, so without compileDirac's nEmptyLoops restoration the loop contributes
   nothing and the kernel comes out 4x too small per collapsed loop. (The dressed path already
   compensated for this; the bare path did not.) Emitted here and graded analytically in
   compare_multiterm_proj_num.cpp against 4 (tr[1]) * 3 (colour Nc) * (q1.ql). *)
netD = ntDeltaDirac[m1, m2] ntDeltaDirac[m2, m1] *
       ntSUNDeltaFund[3, n1, n2] ntSUNDeltaFund[3, n2, n1] *
       ntVec[Global`q1, r1] ntVec[Global`ql, r1];
ntkD = NumTrace[netD, "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
If[Length[ntkD[[1]]["Diagrams"]] < 1,
  Print["FAIL: the collapsed-delta-loop diagram did not survive"]; Exit[1]];
MakeNTKernel[ntkD,
  FileNameJoin[{$gen, "gen_deltaloop_num.cpp"}],
  FileNameJoin[{$gen, "Deltaloop_num_kernel.hh"}],
  FileNameJoin[{$gen, "Deltaloop_num_kernels.hh"}],
  "Name" -> "Deltaloop_num_kernel", "Namespace" -> "deltaloop_num", "Dressings" -> {}];
Print["[gen] ok: collapsed delta-loop kernel emitted"];

Print["[gen] Numeric backend ..."];
MakeNTKernel[ntk,
  FileNameJoin[{$gen, "gen_multiterm_proj_num.cpp"}],
  FileNameJoin[{$gen, "Multiterm_proj_num_kernel.hh"}],
  FileNameJoin[{$gen, "Multiterm_proj_num_kernels.hh"}],
  "Name" -> "Multiterm_proj_num_kernel", "Namespace" -> "multiterm_proj_num",
  "Dressings" -> {}];

Print["Multi-term Dirac projector codegen kernels generated."];
