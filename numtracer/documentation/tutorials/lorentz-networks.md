# Lorentz networks

A network of metrics, momentum vectors and transverse projectors has no "components" once you contract it — the numeric engine sums the shared Lorentz indices away over the loop frame and collapses the network to a polynomial in the frame's momentum components, carrying any surviving inverse propagators $1/l^2$ as separate **atoms**. This contracts the free gamma legs from the [previous tutorial](dirac-traces.md).
We apply it to the warm-up $p\cdot P(l)\cdot p$ and recover the angular factor $1-\cos^2\theta$.

## The program

`Tutorials/03-lorentz-networks/lorentz_networks.cpp`:

```cpp
#include <numtracer.hpp> // the whole NumTracer API — here: nvec / nprojT (NNet) + numeric_value -> MPoly

#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

namespace nm = numtracer::numeric;
namespace net = numtracer::network;
using numtracer::Cx;

// Name the Lorentz axis labels with a single unscoped enum (every label auto-numbered, so all
// distinct); the momentum ids (p = 0, l = 1, indexing the component table) stay plain integers.
enum { mu, nu }; // Lorentz indices

int main() {
  // The loop frame: each momentum is a vid into a component table comp[vid][0..3], whose entries are
  // MPolys. We pick the concrete one-angle frame up front — p along axis 0, l in the 0-1 plane — and
  // keep only its NON-ZERO components symbolic (one scalar variable each), so numeric_value contracts
  // the network into the polynomial in exactly those frame scalars:
  //   p_0 = var 0,   l_0 = var 1,   l_1 = var 2.
  const int nsym = 3;
  nm::LorentzEnv env(nsym); // bind the 3-symbol space once; every MPoly below is minted from it
  std::vector<std::array<nm::MPoly, 4>> comp(2, {env.zero(), env.zero(), env.zero(), env.zero()});
  comp[0][0] = env.var(0); // p = (p_0, 0, 0, 0)
  comp[1][0] = env.var(1); // l = (l_0, l_1, 0, 0)
  comp[1][1] = env.var(2);

  // The transverse projector P(l)_{mu nu} = delta_{mu nu} - l_mu l_nu / l^2 carries its 1/l^2 as
  // "atom" id 0; numeric_value needs that atom's denominator l^2 = sum_i comp[1][i]^2.
  nm::MPoly l2 = env.zero();
  for (int i = 0; i < 4; ++i) l2 = l2 + comp[1][i] * comp[1][i];
  std::vector<nm::MPoly> atomDen = {l2};

  // The Lorentz network as one product term (coefficient 1) of three factors. Sharing index mu
  // (p with the projector) and index nu (projector with the second p) sums both Lorentz indices
  // away -> a scalar network: p.P(l).p.
  nm::NNet lor = {nm::NTerm{Cx{1, 0}, {nm::nvec(mu, {{1.0, 0}}),
                                       nm::nprojT(mu, nu, {{1.0, 1}}, 0),
                                       nm::nvec(nu, {{1.0, 0}})}}};

  // Contract: empty Dirac chain (this is a pure-Lorentz network). The result is one MPoly =
  // p.P(l).p = sp(p,p) - sp(p,l)^2 / l^2 = p_0^2 - p_0^2 l_0^2 / l^2 — two monomials in this frame
  // (the second carries the 1/l^2 atom).
  nm::MPoly poly = env.numeric_value(net::DiracNet{}, lor, comp, atomDen);

  // Evaluate the contracted polynomial at concrete values: p along axis 0, l at angle theta.
  const double Pm = 1.3, l0 = 0.5, l1 = 0.7;
  const std::vector<double> x = {Pm, l0, l1}; // var 0 = p_0, var 1 = l_0, var 2 = l_1
  const double l2v = l0 * l0 + l1 * l1;
  const std::vector<double> atomVal = {1.0 / l2v}; // value of atom 0 = 1/l^2 at this frame
  const double val = nm::eval(poly, x, atomVal).re;

  const double cth = l0 / std::sqrt(l2v); // cos theta
  std::printf("contracted monomials = %zu   (p.P.p = sp(p,p) - sp(p,l)^2 / l^2)\n", poly.size());
  std::printf("p.P(l).p             = %g   (= p^2 (1 - cos^2) = %g)\n", val, Pm * Pm * (1 - cth * cth));
  std::printf("p.P(l).p / p^2       = %g   (= 1 - cos^2 theta = %g)\n", val / (Pm * Pm), 1 - cth * cth);

  const bool ok = std::fabs(val - Pm * Pm * (1 - cth * cth)) < 1e-12;
  std::printf(ok ? "ALL TESTS PASSED\n" : "TESTS FAILED\n");
  return ok ? 0 : 1;
}
```

```bash
cmake --build build --target lorentz_networks && ./build/lorentz_networks
```
```text
contracted monomials = 2   (p.P.p = sp(p,p) - sp(p,l)^2 / l^2)
p.P(l).p             = 1.11905   (= p^2 (1 - cos^2) = 1.11905)
p.P(l).p / p^2       = 0.662162   (= 1 - cos^2 theta = 0.662162)
ALL TESTS PASSED
```

That `0.662162` is the angular factor $1-\cos^2\theta$ for this frame — reached by contracting the
network numerically over the loop frame, with no component-by-component bookkeeping.

## Reading it

An `NNet` is a Lorentz network: a sum of product terms (`NTerm`), each a coefficient times a list
of factors. The factor builders are:

| builder                       | meaning                                                                                                                                                |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `nvec(Lbl, vlc)`              | a 4-vector on Lorentz index `Lbl`, momentum $=\sum \text{coeff}\cdot\text{comp}(vid)$ (`vlc` = list of `{coeff, vid}`)                                  |
| `nmet(Mu, Nu)`                | the metric $\delta_{\mu\nu}$                                                                                                                            |
| `nprojT(Mu, Nu, vlc, atom)`    | the transverse projector $P^T_{\mu\nu}(l)=\delta_{\mu\nu}-l_\mu l_\nu/l^2$; `l` $=\sum\text{coeff}\cdot\text{comp}(vid)$, `atom` is the id of its $1/l^2$ |
| `nprojL(Mu, Nu, vlc, atom)`    | the longitudinal projector $P^L_{\mu\nu}(l)=l_\mu l_\nu/l^2$ (same arguments as `nprojT`); note $P^T+P^L=\delta$                                        |
| `neps(A, B, C, D)`            | the Levi-Civita tensor $\varepsilon_{ABCD}$ (the $\gamma_5$ trace's antisymmetric tensor)                                                              |

`env.numeric_value(dirac, net, comp, atomDen)` — a method on the `LorentzEnv` that carries `nsym` —
contracts the Dirac chain (empty here) against the
Lorentz network and returns one `MPoly` — a polynomial in the frame's component variables, each
monomial optionally carrying inverse-propagator **atoms**. For $p\cdot P(l)\cdot p$ it has **two**
monomials: $\mathrm{sp}(p,p)$ and $-\mathrm{sp}(p,l)^2\cdot(1/l^2)$ — exactly
$p^2 - (p\cdot l)^2/l^2 = p^2(1-\cos^2\theta)$. `eval(poly, x, atomVal)` then substitutes the
component values `x` and the atom values `atomVal` to get the number. The contraction stays bounded
— it never expands the projector mask into a $2^k$ component blow-up — which is what keeps the
kernel small and fast. See [the numeric engine](../internals/numeric-engine.md).

Next, [a full diagram](full-diagram.md): the projector contracting the Dirac trace's free legs,
assembled and validated end to end.
