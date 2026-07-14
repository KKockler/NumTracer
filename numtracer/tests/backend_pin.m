(* ::Package:: *)

(* backend_pin.m — TEST fixture. Optional override of FunKit's evaluation backend.

   By DEFAULT this does nothing: the tests run on whatever backend FunKit selects, which since
   b520d21 is the C++ one. That is deliberate — it is the backend production flows use, so it is the
   one the suite should exercise.

   Why the file exists at all. FunKit's `FRoute` picks the loop momentum from the leg order it is
   handed, and the two backends order identical-field vertices differently. So the backend choice used
   to change the loop-momentum ROUTING, and hence the integrand at fixed |l| — while leaving the
   integral (the physics) untouched. Any test that graded on the POINTWISE integrand therefore failed
   by O(1) on a backend switch with nothing actually wrong. Two things fixed that:

     * FunKit's routing algorithm was reworked (2026-07-14) so both backends now route alike (up to
       the sign convention of the loop momentum);
     * the compare harnesses grade on the ANGULAR-INTEGRATED value, which is invariant under the
       remaining reparametrisation freedom (see compare_zaqbq1.cpp / compare_za4_num.cpp).

   With both in place the pin is no longer needed, so it is off. It is kept as an escape hatch: if a
   future FunKit change makes the routing diverge again and you need to bisect against the frozen
   oracles, force the old behaviour with

       NT_FUNKIT_BACKEND=Mathematica   (or =Cpp)

   Reminder: if you ever re-freeze the test artifacts, regenerate BOTH the numeric kernels and the
   refshim FORM oracles from the same run — see tests/gen/README.md. *)

With[{nt$be = Environment["NT_FUNKIT_BACKEND"]},
  Which[
    nt$be === "Mathematica", FSetBackendMathematica[],
    nt$be === "Cpp",         FSetBackendCpp[],
    True,                    Null  (* default: leave FunKit's own choice alone *)
  ]
];
