#!/usr/bin/env wolframscript
(* Codegen gate for PROJECTOR ORTHONORMALITY through the full NumTracer pipeline.
   For a tensor basis, the dual projectors must satisfy  <P_i, T_j> = KroneckerDelta[i,j]  exactly and
   independently of the kinematics. This is the guard that would have caught the multi-term-projector
   bug (an eager Dirac Plus of mixed gamma parity traced to an identically-zero kernel) at the level
   of real physics bases rather than a synthetic net.

   WEIGHTED TRACE. One kernel per pair would mean n^2 committed headers per basis (64 for
   AqbqDirect8). Instead contract the WHOLE matrix against a fixed weight matrix in a single net:
       W = Sum_ij w[i,j] <P_i, T_j>     must equal     Sum_i w[i,i]
   The weights are distinct primes, so no plausible error cancels: any single wrong entry shifts W by
   its own weight, and the diagonal/off-diagonal weights are disjoint. This covers every one of the
   n^2 pairs in ONE kernel. compare_ortho_num.cpp grades W against Sum_i w[i,i] at random momenta;
   the per-pair matrix (for diagnosing WHICH entry moved) is the opt-in sweep, not this gate.

   CONTRACTION. The dual projector of a VERTEX basis is dual to the VERTEX, so the relation to test
   is <P_i, V_j> = delta_ij with V_j = TBGetVertex. Do NOT route this through TBGetInnerProduct:
   TBMakeInnerProduct selects TBCanonicalProduct for a vertex basis (AqbqDirect is one), i.e. the
   product used to project onto vertices, so pairing it with TBGetBasisElement mixes conventions.
   Contract directly, with the indices supplied by hand.

   INDEX CONVENTION (legs of AqbqDirect: 1 = gluon {p1,mu,a1}, 2 = anti-quark {p2,d2,A2,F2},
   3 = quark {p3,d3,A3,F3}). The quark leg of the projector closes against the ANTI-quark leg of the
   vertex, and the MOMENTUM TRAVELS WITH ITS LEG — the whole bundle moves. Swapping only the
   (Dirac,colour,flavour) triples while pinning p2/p3 transposes the Dirac structure instead, which
   shows up as sign flips on the C-odd structures and a 2<->3 leak. *)
$repo = AbsoluteFileName @ FileNameJoin[{DirectoryName[$InputFileName], "..", "..", ".."}];
$gen  = FileNameJoin[{$repo, "numtracer", "tests", "gen"}];
If[Quiet[Needs["FunKit`"]] === $Failed,
  $Path = Prepend[$Path, FileNameJoin[{$HomeDirectory, "Documents", "Uni", "Code", "FunKit"}]]];
Needs["FunKit`"];
Get["DiFfRG`"];
FSetRegisterSize[64];
SetDirectory[$repo];   (* qcd147_setup.m resolves its TBCache relative to the repo root *)
Get[FileNameJoin[{$repo, "numtracer", "tests", "qcd147_setup.m"}]];
Get[FileNameJoin[{$repo, "numtracer", "mathematica", "NumTracer.m"}]];
SetNc[3]; SetNf[2];


SPParam[expr_] := UseLorentzLinearity[expr] //. {
   lf1 -> l1, sp[p, p] -> p^2, sp[l1, l1] -> l1^2,
   sp[l1, p1] -> p*l1*cos[l1, p1], sp[l1, p2] -> p*l1*cos[l1, p2],
   sp[l1, p3] -> p*l1*cos[l1, p3],
   Sqrt[(a_)^2] :> a, ((a_)^2)^((n_)/2) :> a^n,
   Power[Power[l1_, 2], Rational[n_, 2]] :> l1^n,
   cos[l1, p1] :> cosl1p1, cos[l1, p2] :> cosl1p2, cos[l1, p3] :> cosl1p3};

(* GENERAL 3-point kinematics, NOT the symmetric point. sp3Frame pins a measure-zero slice on which
   the full AqbqDirect basis is degenerate (Det[Gram] = 0 there, so every dual projector is 0/0);
   checking orthonormality only at that point would test the tracer where the algebra breaks down.
   Verified: Det[Gram] = 0 at (1,1,-1/2) and 1.87e22 at (1,7/4,-3/10). *)
frame = NumTracer`gen3Frame[p1m, p2m, cosP, l1, cos1, cos2, p1, p2, p3, l1];
args  = {l1, cos1, cos2, p1m, p2m, cosP};

(* distinct primes as weights: w[i,j] = prime #(n*(i-1)+j) *)
wOf[n_][i_, j_] := Prime[n*(i - 1) + j];

idxP = {{p1, mu, a1}, {p2, d2, A2, F2}, {p3, d3, A3, F3}};
idxV = {{p1, mu, a1}, {p3, d3, A3, F3}, {p2, d2, A2, F2}};

doBasis[b_String, nm_String, ns_String, n_Integer] := Module[{net, ntk, expect, w},
  (* Guard against a VACUOUS pass: if the basis never resolved (no TBCache here and FORM did not
     build it), TBGetProjector stays unevaluated, the net degenerates to a symbolic expression, and
     NumTrace happily returns a meaningless diagram. Assert the accessors actually evaluated. *)
  If[! FreeQ[TBGetProjector[b, 1, idxP], _TBGetProjector] ||
     ! FreeQ[TBGetVertex[b, 1, idxV], _TBGetVertex],
    Print["FAIL: basis ", b, " did not resolve (TBGetProjector/TBGetVertex unevaluated). ",
          "Without a TBCache in tests/gen, FORM must build the Gram — check the basis setup ",
          "rather than trusting this gate."]; Exit[1]];

  expect = Total[Table[wOf[n][i, i], {i, n}]];
  Print["[ortho] ", b, ": n=", n, "  expected weighted trace = ", expect];
  w   = Table[wOf[n][i, j], {i, n}, {j, n}];
  net = Sum[w[[i, j]] * FTerm[TBGetProjector[b, i, idxP], TBGetVertex[b, j, idxV]], {i, n}, {j, n}];
  ntk = NumTrace[FromFunKit[SPParam[net]], "Frame" -> frame, "Args" -> args, "Dressings" -> {}];
  If[Length[ntk[[1]]["Diagrams"]] < 1,
    Print["FAIL: ", b, " weighted orthonormality trace produced NO diagrams — an identically-zero ",
          "kernel. This is the silent-zero signature (a multi-term projector dropped by the ",
          "odd-trace filter, or a component misrouted to the colour compiler)."]; Exit[1]];
  Print["[ortho] ", b, ": ", Length[ntk[[1]]["Diagrams"]], " diagram(s)"];
  MakeNTKernel[ntk,
    FileNameJoin[{$gen, "gen_" <> ns <> ".cpp"}],
    FileNameJoin[{$gen, nm <> ".hh"}],
    FileNameJoin[{$gen, nm <> "s.hh"}],
    "Name" -> nm, "Namespace" -> ns, "Dressings" -> {}];
  Print["[ortho] ", b, ": kernel ", nm, " emitted (expect ", expect, ")"];
  expect];

e147 = doBasis["AqbqDirect147", "Ortho147_kernel", "ortho147", 3];
(* The FULL basis: no restriction, no TBCache — it ships with TensorBases. 12 structures, an oblique
   Gram (16 nonzero off-diagonals), all 144 pairs in one weighted net. Only tractable in the GENERAL
   frame: at the symmetric point its Gram is singular and every projector degenerates to 0/0. *)
e12  = doBasis["AqbqDirect", "Ortho12_kernel", "ortho12", 12];

Print["Orthonormality codegen kernels generated. AqbqDirect147 -> ", e147, ", AqbqDirect(full) -> ", e12];
