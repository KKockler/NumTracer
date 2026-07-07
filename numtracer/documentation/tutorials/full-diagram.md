# A full diagram

The previous tutorials built each piece in isolation: the [colour factor](color-contraction.md),
the [Dirac trace](dirac-traces.md), the [Lorentz network](lorentz-networks.md). This capstone
assembles the whole quark-self-energy numerator,

$$
T_{\text{num}} = \mathrm{tr}\!\big[\slashed p\,\gamma^\mu\,\slashed q\,\gamma^\nu\big]\,
P_{\mu\nu}(l), \qquad q = l - p,
$$

**by hand** with the numeric engine, and validates it two ways at one frame point: the numeric
contraction and the closed form.

## The program

`Tutorials/04-full-diagram/full_diagram.cpp` (the numeric path; the closed-form check follows in
the source):

```cpp
#include "numtracer/numeric/numeric_contract.hpp" // numeric_value, nprojT, dslash, dgamma
// … <array> <cmath> <cstdio> <vector> …

namespace nm = numtracer::numeric;
namespace net = numtracer::network;
using numtracer::Cx;

// One-angle frame: p along axis 0, l at angle theta in the 0-1 plane, q = l - p.
const double p = 1.3, l0 = 0.5, l1y = 0.7;
const double pvec[4] = {p, 0, 0, 0}, lvec[4] = {l0, l1y, 0, 0}, qvec[4] = {l0 - p, l1y, 0, 0};
const double l2 = l0 * l0 + l1y * l1y;

enum { mu, nu }; // name the free gamma legs' Lorentz indices (one unscoped enum, all distinct)

// The numeric path. Momenta are numbers here, so the components are CONSTANT polynomials
// (nsym = 0 symbols). comp[vid][m] = component m of momentum vid: p -> 0, q -> 1, l -> 2.
const int nsym = 0;
std::vector<std::array<nm::MPoly, 4>> comp(3);
// … fill comp[0]=p, comp[1]=q, comp[2]=l with MPoly::constant …

// The closed chain p/ gamma^mu q/ gamma^nu, the two free gammas on Lorentz legs mu, nu.
net::DiracNet chain = {net::dslash({{1.0, 0}}), net::dgamma(mu),
                       net::dslash({{1.0, 1}}), net::dgamma(nu)};
// Lorentz net = the transverse projector P_{mu nu}(l): legs mu/nu, momentum vid 2, and its
// 1/l^2 factor tracked as inverse-atom id 0 (value supplied in atomDen).
nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nprojT(mu, nu, {{1.0, 2}}, 0)}}};
std::vector<nm::MPoly> atomDen = {/* l^2 as an MPoly */};

nm::MPoly tr = nm::numeric_value(nsym, chain, lor, comp, atomDen);
Cx num = nm::eval(tr, /*symbols*/ {}, /*atom values*/ {1.0 / l2}); // -> T_num
```

```bash
cmake --build build --target full_diagram && ./build/full_diagram
```
```text
numeric  T_num = 3.527568
closed   T_num = 3.527568   (4 p(-3 c l1 + p + 2 c^2 p), q = l - p)
ALL TESTS PASSED
```

## Reading it

The diagram is one Dirac chain contracted against one Lorentz network — exactly the pieces from
tutorials 2 and 3, now joined:

| piece | meaning |
|---|---|
| `comp[vid][mu]` | component $\mu$ of momentum `vid`; here constant `MPoly`s (numbers), so `nsym = 0` |
| `dslash({{1.0, 0}})`, `dslash({{1.0, 1}})` | the slashed momenta $\slashed p$, $\slashed q$ |
| `dgamma(mu)`, `dgamma(nu)` | the free gammas $\gamma^\mu$, $\gamma^\nu$, open Lorentz legs `mu`/`nu` |
| `nprojT(mu, nu, {{1.0,2}}, 0)` | the projector $P_{\mu\nu}(l)$: legs `mu`/`nu`, momentum `vid 2` ($l$), inverse-atom id 0 |
| `atomDen = {l²}` | tells the engine atom 0 is $l^2$, so `eval` can supply $1/l^2$ |
| `numeric_value(...)` | folds the Dirac chain (matrix products) against the Lorentz net (index elimination) → one polynomial |

The momentum components are plain numbers here, so `numeric_value` returns a constant `MPoly` and
`eval` reads off the value. (In a generated kernel they are instead polynomials in the frame
variables $|l|, \cos\theta, |p|$, and `eval` is replaced by the lowered straight-line code — same
contraction, symbolic inputs.)

Two independent computations agreeing to machine precision —
$\text{numeric} = \text{closed} = 3.527568$ — is what pins the algebra down: the closed form
$4p(-3c\,l_1 + p + 2c^2 p)$ is the analytic value for $q = l - p$, and the numeric engine reaches
the same number the way the code generator does.

## Lowering to an optimized kernel

The contraction above is only the first half of the pipeline. The second half — **lowering** the
resulting polynomial to straight-line real arithmetic and **emitting** it as C++ — is callable right
here from the same program, no Mathematica involved. This is a preview of the codegen; the passes
below (Horner factoring, real CSE) are explained in full in
[CSE and Horner lowering](../internals/cse-and-lowering.md). It produces the *same bytes* the
[front-end](generating-kernels.md) writes:

```cpp
#include "numtracer/codegen/gen.hpp" // GlobalEnv, GenProg, emit_cpp, emit_env_layout

net::GlobalEnv g;                          // the shared symbol table (one f[] slot per symbol/atom)
net::GenProg  prog = nm::to_genprog(tr, g); // Horner-factor + hash-consed real CSE
net::emit_env_layout(std::cout, g);         // comment: which f[i] is which symbol / 1/l^2
net::emit_cpp(std::cout, prog, "T_num");    // -> double T_num(const double* f) { ... }
```

One catch: in the contraction above the momenta are plain numbers, so `tr` is a **constant** `MPoly`
and would lower to a bare number. To get a real *kernel* the tutorial redoes the contraction with the
momentum **components as symbols** — `ns = 3`: `f[0] = |p|`, `f[1] = l_0`, `f[2] = l_{1y}` — reusing
the *same* `chain`/`lor` tokens (they refer to momentum ids, not values). The emitted function is then
a genuine function of the frame:

```text
// fundamental-symbol env layout (fill f[i] per call):
//   f[0] = var(0)   // |p|
//   f[1] = var(1)   // l0
//   f[2] = var(2)   // l1y
//   f[3] = inv(0)   // 1 / l^2
static inline double T_num(const double *f) {
  const double s0 = -4;
  const double s1 = 4;
  const double s2 = f[0];
  const double s3 = s1*s2;
  // … Horner-factored, every repeated subexpression shared exactly once …
  const double s23 = s2*s22;
  return s23;
}
```

That is exactly the `trN(const double* f)` shape of a generated kernel: a flat straight-line body over
a shared `f[]` environment, which a companion `fill(...)` populates from the kernel's scalar
arguments once per call.

**Can you run the lowered kernel directly?** The lowered form *is* C++ — you run it by compiling
`T_num(f)` and calling it (that is what a generated kernel does; there is no interpreter for the
lowered instruction stream). In-process, what you can run without a compile step is the *contracted
polynomial* itself: `nm::eval(tr, symbolValues, atomValues)` evaluates it directly — the same numbers,
just not the Horner/CSE'd arithmetic. The tutorial uses exactly that to self-check the symbolic
contraction against the closed form before printing the kernel.

## From here to a kernel

A real flow has many such diagrams, each multiplied by its scalar coefficient (dressings,
regulators, propagator denominators), summed. Writing all of that by hand is mechanical and
derivable from the flow equation — which is what the [Mathematica front-end](generating-kernels.md)
automates, emitting a flat C++ kernel from a tensor network in a small DSL.
