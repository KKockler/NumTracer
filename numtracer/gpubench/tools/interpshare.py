#!/usr/bin/env python3
"""interpshare.py — share the interpolator coordinate transform across dressings.

The measured cost of one SplineInterpolator1D lookup in these kernels is ~210 fp64 instructions,
and ~200 of those are a single fp64 log1p inside LogarithmicCoordinates1D::backward. A generated
kernel evaluates many dressings at the SAME momentum — nf2 ZAqbq4 does 347 lookups over only 145
distinct arguments — and pays that log1p once per lookup, because each interpolator owns its own
`coordinates` members and nvcc cannot prove they are equal across objects, so it cannot CSE.

This rewrites the emitted wrapper so each distinct argument computes the grid index once:

    const auto _interp40 = ZAqbq1(<arg>);          ->   const auto _ix3 = ZAqbq1.index(<arg>);
    const auto _interp48 = ZAqbq4(<arg>);               const auto _interp40 = ZAqbq1.at(_ix3);
    const auto _interp50 = ZAqbq7(<arg>);               const auto _interp48 = ZAqbq4.at(_ix3);
                                                        const auto _interp50 = ZAqbq7.at(_ix3);

CORRECTNESS PRECONDITION: the index is only transferable between interpolators that share a
coordinate system. In these models every dressing spline is constructed on the same
`coordinates1D` (model.hh builds them all from one Coordinates1D), so this holds — but it is an
assumption about the consumer, not a property the rewrite can check. The index is taken from the
first interpolator that uses the argument, and the correctness gate (checksum vs baseline) is what
actually confirms it.

Regulator calls (RB/RBdot/RF/RFdot/...) are left alone: they are polynomials, not table lookups.

Usage:  interpshare.py <in_kernel.hh> -o <out_kernel.hh> [--report]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

RE_INTERP = re.compile(r"^(\s*)const auto (_interp\d+) = (\w+)\((.*)\);\s*$")
# The kernel signature is the ground truth for which callees are actually splines. Guessing from a
# denylist is not enough: the same `const auto _interpN = f(...)` form also carries regulator calls
# (RB/RFdot/...) and, on complex flows, ntRe(ns::trK(fenv)) — and `ntRe` is an overloaded member
# function, so rewriting it to `ntRe.index(...)` does not even compile.
# NB: the template argument list is nested (SplineInterpolator1D<double,
# LogarithmicCoordinates1D<double>, GPU_memory>&), so a [^>]* body stops at the inner '>'. Scan to
# the first '&' instead — a reference marker only appears at the end of each parameter's type.
RE_SPLINE_PARAM = re.compile(
    r"(?:SplineInterpolator1D|LinearInterpolator1D)\s*<[^&]*?&\s*(\w+)")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("header", type=Path)
    ap.add_argument("-o", "--out", type=Path)
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()

    text = args.header.read_text(errors="replace")
    lines = text.splitlines()

    # The kernel's own signature names every interpolator it was handed. Anything else appearing in
    # the `const auto _interpN = f(...)` form is a regulator or ntRe/ntIm and must be left alone.
    splines = set(RE_SPLINE_PARAM.findall(text))
    if not splines:
        print(f"{args.header.parent.name:<12} no interpolator parameters — nothing to do")
        if args.out and not args.report:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(text)
        return 0

    # first pass: which arguments repeat across interpolator lookups?
    seen: dict[str, int] = {}
    for ln in lines:
        m = RE_INTERP.match(ln)
        if not m or m.group(3) not in splines:
            continue
        seen[m.group(4).strip()] = seen.get(m.group(4).strip(), 0) + 1

    shared = {a for a, n in seen.items() if n >= 2}
    idx_name = {a: f"_ix{i}" for i, a in enumerate(sorted(shared))}
    emitted: set[str] = set()

    out: list[str] = []
    n_split = 0
    for ln in lines:
        # Complex flows emit three alternative bodies under #if <VERDICT>==2 / ==1 / #else. Those
        # branches are mutually exclusive, so an index hoisted in one is not in scope in another —
        # reset the emitted set at every preprocessor boundary and let each branch own its hoists.
        if ln.lstrip().startswith(("#if", "#else", "#elif", "#endif")):
            emitted.clear()
            out.append(ln)
            continue
        m = RE_INTERP.match(ln)
        if not m or m.group(3) not in splines:
            out.append(ln)
            continue
        indent, name, callee, arg = m.group(1), m.group(2), m.group(3), m.group(4).strip()
        if arg not in shared:
            out.append(ln)
            continue
        ix = idx_name[arg]
        if arg not in emitted:
            # the first user of this argument owns the transform; every later one reuses the index
            out.append(f"{indent}const auto {ix} = {callee}.index({arg});")
            emitted.add(arg)
        out.append(f"{indent}const auto {name} = {callee}.at({ix});")
        n_split += 1

    total = sum(seen.values())
    print(f"{args.header.parent.name:<12} splines={total:<4} distinct={len(seen):<4} "
          f"shared_args={len(shared):<4} rewritten={n_split:<4} "
          f"log1p {total} -> {len(seen)} ({100 * (total - len(seen)) // max(1, total)}% saved)")

    if args.out and not args.report:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
