#!/usr/bin/env wolframscript
(* step-10 — From the Wetterich equation to a NumTracer network.
 *
 * Steps 06-09 handed NumTrace a network we wrote by hand. In a real fRG calculation nobody writes
 * the network: it is DERIVED from the flow equation. This script does the derivation for the
 * running example — the quark two-point flow, whose numerator steps 03-05 built by hand — and hands
 * the result straight to NumTrace.
 *
 * TIER C. Needs the full fRG toolchain: FunKit (the functional-derivative engine), TensorBases (the
 * tensor bases and projectors), DiFfRG, and FORM. Unlike steps 06-09 this is NOT wired into ctest —
 * it is here to be read and, if you have the toolchain, run:
 *
 *     NT_TUTORIAL_OUT=/tmp/out wolframscript -script quark_selfenergy.wl
 *
 * It writes its kernel to NT_TUTORIAL_OUT and touches nothing in the repository.
 *)

(* --------------------------------------------------------------------------------------------
 * @snip begin: load
 * The toolchain. FunKit supplies WetterichEquation / FTakeDerivatives / FTruncate / FRoute and the
 * COEN C++ emitter; DiFfRG supplies the regulator vocabulary; TensorBases (loaded by FunKit)
 * supplies TBGetProjector. FORM is used by TensorBases to compute the tensor BASES — never to do
 * NumTracer's trace.
 *)
If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"]; Get["DiFfRG`"]; FSetRegisterSize[64];

$repo = AbsoluteFileName@FileNameJoin[{DirectoryName[$InputFileName], "..", ".."}];
SetDirectory[$repo];                                   (* TBCache/ is resolved relative to here *)
DefineFormExecutable["/usr/bin/tform -w16"];

(* The per-theory setup: field content, truncation, tensor bases and dressing rules. This is NOT
   part of NumTracer — it is what YOUR model file provides. numtracer/tests/qcd_setup.m is a
   trimmed copy of the production QCD_Nf2 setup; read it alongside this script. *)
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
(* @snip end: load *)

$out = Environment["NT_TUTORIAL_OUT"];
If[$out === $Failed, $out = Directory[]];

(* --------------------------------------------------------------------------------------------
 * @snip begin: derive
 * (1) DIFFERENTIATE. The Wetterich equation is an identity for the scale derivative of the
 * effective action. Two functional derivatives with respect to qb and q give the flow of the quark
 * TWO-POINT function; the external index labels i1, i2 are carried along symbolically.
 *
 *   FTruncate  drops everything outside the truncation declared in the setup.
 *   FRoute     assigns loop-momentum routing to each diagram — this is where "which line carries
 *              l, which carries l-p" is decided, and it is leg-order sensitive (see backend_pin.m).
 *)
fRGqbq = FTakeDerivatives[WetterichEquation, {qb[i1], q[i2]}] // FTruncate // FRoute;

Print["step-10: 1-loop external indices = ", fRGqbq["1-Loop"]["ExternalIndices"]];
(* @snip end: derive *)

(* --------------------------------------------------------------------------------------------
 * @snip begin: project
 * (2) PROJECT. The flow of a two-point function is still a matrix in Dirac and colour space. To get
 * a SCALAR flow we contract it with the dual of the basis element we want — here element 1 of the
 * "qbq" basis, the quark wave function.
 *
 * Three details, each of which has cost somebody a week:
 *   - `**` is FunKit's NON-COMMUTATIVE product. The projector must multiply the loop from the left.
 *   - The index order is {i2, i1} — REVERSED with respect to the derivative order {qb[i1], q[i2]}.
 *   - The leading MINUS compensates the sign convention of the qbq dressing (see qcd_setup.m,
 *     where dressing[GammaN, {qb, q}, 1, ...] :> -Zq[...]). Drop it and the whole flow flips sign.
 * step-11 is about getting these right, and about asserting them instead of hoping.
 *)
trace = FTerm[-TBGetProjector["qbq", 1, {i2, i1} /. fRGqbq["1-Loop"]["ExternalIndices"]]] **
        (fRGqbq["1-Loop"]["Expression"] /. FMakeDiagrammaticRules[]);
(* @snip end: project *)

(* --------------------------------------------------------------------------------------------
 * @snip begin: dress
 * (3) DRESS and REDUCE. Two substitution layers, in this order:
 *
 *   dressingRules  maps the abstract `dressing[GammaN, {A,A}, 1, {...}]` heads produced by the
 *                  derivative onto the NAMED functions your kernel will take as parameters
 *                  (ZA, Zq, Mq, ...), and inserts the regulator.
 *   PropParam      does the SCALAR kinematic reduction: sp[l1,p1] -> l1 p cos1, etc. This is the
 *                  scalar half only — the tensor half gets its kinematics from the frame.
 *
 * A `dressing[...]` head that no rule matches does NOT fail: it survives into the integrand and is
 * emitted into the kernel as an undefined `dressing(GammaN, List(...))` call. See step-22.
 *)
flow = PropParam[trace // dressingRules];
(* @snip end: dress *)

(* --------------------------------------------------------------------------------------------
 * @snip begin: import
 * (4) IMPORT. FromFunKit rewrites the traced FunKit flow into the NumTracer DSL: gammas become
 * ntGamma, momentum legs become ntVec, scalar products become ntSP, and SU(N) tokens get their rank
 * baked in (colour from Nc, isospin from "FlavourGroup"). From here on it is exactly the object
 * steps 06-09 built by hand.
 *)
net = FromFunKit[flow];
(* @snip end: import *)

(* --------------------------------------------------------------------------------------------
 * @snip begin: emit
 * (5) TRACE and EMIT — identical to step-06.
 *)
frame = NumTracer`propFrame[Global`p, Global`l1, Global`cos1, Global`p1, Global`l1];
args  = {Global`l1, Global`cos1, Global`p, Global`k};
dress = {ZA, Zq, Mq, ZAqbq1, dtZA, dtZq};

ntk = NumTrace[net, "Frame" -> frame, "Args" -> args];
Print["step-10: diagrams = ", Length[ntk[[1]]["Diagrams"]], "   (expected 2: the regulator ",
      "insertion sits on the gluon line in one and on the quark line in the other)"];
If[Length[ntk[[1]]["Diagrams"]] =!= 2,
  Print["step-10 FAIL: expected 2 diagrams"]; Exit[1]];

MakeNTKernel[ntk,
  FileNameJoin[{$out, "gen_step10.cpp"}],
  FileNameJoin[{$out, "Zq_tutorial_kernel.hh"}],
  FileNameJoin[{$out, "Zq_tutorial_kernels.hh"}],
  "Name" -> "Zq_tutorial_kernel", "Namespace" -> "step10", "Dressings" -> dress];
(* @snip end: emit *)

Print["step-10: wrote Zq_tutorial_kernel.hh to ", $out];
