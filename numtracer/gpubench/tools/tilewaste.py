#!/usr/bin/env python3
"""tilewaste.py — how many of the threads an MDRange launch starts actually do work.

Kokkos' MDRangePolicy tiles the iteration space and launches one block per tile, so a tile
dimension that does not DIVIDE its extent launches masked-off lanes. `Kokkos::LaunchBounds<maxT>`
clamps the tile product, which means DIFFRG_LAUNCH_BOUNDS silently sets the block SHAPE as well as
hinting registers to ptxas -- and the A100 LaunchBounds sweep therefore measured the sum of two
unrelated effects while attributing it to occupancy alone.

This script computes the masked fraction per flow, statically. The tile is chosen host-side from
the extents and maxTperB only, with no GPU-model dependence, so the answer is exact and needs no
A100 -- it is the same number the cluster would print.

The tiling logic below is a transcription of Kokkos' own
`MDRangePolicy::update_tiling_properties` (bundled/include/KokkosExp_MDRangePolicy.hpp:383-450)
plus `TileSizeRecommended<Cuda>` (bundled/include/Cuda/Kokkos_Cuda_MDRangePolicy.hpp:25-70).
Keep it in sync if the bundled Kokkos is bumped -- `--selftest` checks the parts that are
independently known.

Usage:
  tilewaste.py <app_dir> [<app_dir> ...]      # app dir holds parameter.{json,toml} + flows/
  tilewaste.py --lb 64,96,128,192,256 <dir>   # sweep launch-bounds values
  tilewaste.py --selftest                     # verify against hand-derived cases
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

CUDA_MAX_TOTAL_TILE_SIZE = 512  # get_tile_size_properties<Cuda>: properties.max_total_tile_size

# TileSizeRecommended<Cuda>, Iterate::Left (default_inner_direction<Cuda> == Iterate::Left).
RECOMMENDED_LEFT = {
    2: [64, 4],
    3: [32, 2, 4],
    4: [16, 2, 2, 4],
    5: [16, 2, 4, 2, 1],
    6: [8, 4, 2, 2, 2, 1],
}

# QuadratureIntegrator<dim, ...> and the config keys each integrator reads for its grid, from
# DiFfRG/physics/integration/vacuum/integrator_p2_*.hh. The launched rank is 1 + dim: extents[0] is
# the external grid (integral_view.size()), extents[1..] are the quadrature orders.
INTEGRATORS = {
    "Integrator_p2_4D_3ang": ["x_order", "cos1_order", "cos2_order", "phi_order"],
    "Integrator_p2_4D_2ang": ["x_order", "cos1_order", "cos2_order"],
    "Integrator_p2_1ang": ["x_order", "cos1_order"],
    # Finite-T: QuadratureIntegrator_fT sets grid_size[dim-1] = matsubara_nodes.size() at runtime
    # (quadrature_integrator_fT.hh:63), so the LAST axis is Matsubara, not a config key. It varies
    # with T and matsubara_precision_factor, hence the MATSUBARA sentinel and --nmats.
    "Integrator_fT_p2_4D_2ang": ["x_order", "cos1_order", "phi_order", "MATSUBARA"],
    "Integrator_fT_p2_1ang": ["x_order", "cos1_order", "MATSUBARA"],
}
MATSUBARA = "MATSUBARA"


def recommended_tile(rank: int) -> list[int]:
    if rank in RECOMMENDED_LEFT:
        return list(RECOMMENDED_LEFT[rank])
    tile = [2] * rank
    tile[0] = 16
    return tile


def kokkos_tile(extents: list[int], max_t_per_b: int | None) -> list[int]:
    """Reproduce MDRangePolicy::update_tiling_properties for Iterate::Left on Cuda.

    max_t_per_b is LaunchBounds::maxTperB; None or 0 means "no launch bounds", in which case only
    the 512 cap applies.
    """
    rank = len(extents)
    effective = CUDA_MAX_TOTAL_TILE_SIZE
    if max_t_per_b:
        effective = min(effective, max_t_per_b)

    default_tile = recommended_tile(rank)
    tile = [0] * rank
    prod = 1
    # inner_direction == Iterate::Left -> walk dimensions 0 .. rank-1
    for i in range(rank):
        if prod * default_tile[i] <= effective:
            tile[i] = default_tile[i]
        else:
            d = default_tile[i]
            while d > 1 and prod * d > effective:
                d >>= 1
            tile[i] = d if d > 1 else 1
        prod *= tile[i]
    return tile


def divisors_up_to(n: int, cap: int) -> list[int]:
    return [d for d in range(1, min(n, cap) + 1) if n % d == 0]


def divisibility_tile(extents: list[int], budget: int) -> list[int]:
    """The tile this SHOULD use: every dimension divides its extent, so nothing is masked.

    Kokkos' heuristic only ever halves the recommended tile, so it cannot reach a tile like
    {16,4,2,1,1} on extents {64,32,6,6,6} -- 6 is not a power of two and 4 does not divide it.
    That is the whole bug: the shipped LB=128 picks tile 4 against extent 6.

    Search: among all tiles whose dimensions divide their extents and whose product fits the
    budget, take the largest product (biggest block = fewest launches, and the register file is
    allocated per warp so bigger blocks waste less of it). Ties break toward the Kokkos
    recommendation, which encodes the coalescing preference on dim 0.
    """
    rank = len(extents)
    rec = recommended_tile(rank)
    best = None
    # Dimensions are few and divisor sets are small, so an exhaustive walk is fine; cap each
    # dimension's candidate tile at the budget to keep the product bounded.
    cands = [divisors_up_to(e, budget) for e in extents]

    def rec_distance(tile):
        # prefer shapes close to Kokkos' recommendation, weighted toward the innermost dim
        return sum(abs(tile[i] - rec[i]) * (rank - i) for i in range(rank))

    def walk(i, prod, tile):
        nonlocal best
        if i == rank:
            key = (prod, -rec_distance(tile))
            if best is None or key > best[0]:
                best = (key, list(tile))
            return
        for d in cands[i]:
            if prod * d > budget:
                continue
            tile.append(d)
            walk(i + 1, prod * d, tile)
            tile.pop()

    walk(0, 1, [])
    return best[1] if best else [1] * rank


def analyse(extents: list[int], max_t_per_b: int | None, tile: list[int] | None = None) -> dict:
    if tile is None:
        tile = kokkos_tile(extents, max_t_per_b)
    block = 1
    for t in tile:
        block *= t
    tiles_per_dim = [-(-e // t) for e, t in zip(extents, tile)]  # ceil
    nblocks = 1
    for n in tiles_per_dim:
        nblocks *= n
    useful = 1
    for e in extents:
        useful *= e
    launched = nblocks * block
    # which dimensions are responsible for the waste
    bad = [
        (i, extents[i], tile[i])
        for i in range(len(extents))
        if extents[i] % tile[i] != 0
    ]
    return {
        "tile": tile,
        "block": block,
        "grid": tiles_per_dim,
        "blocks": nblocks,
        "useful": useful,
        "launched": launched,
        "waste": launched / useful,
        "bad_dims": bad,
    }


def read_config(app: Path) -> dict:
    """Pull the integration orders out of parameter.json or parameter.toml.

    Only the flat set of *_order keys is needed, so a regex over the TOML avoids a toml dependency
    (python 3.11+ has tomllib, but this must run under the system python too).
    """
    cfg = {}
    pj = app / "parameter.json"
    pt = app / "parameter.toml"
    if pj.exists():
        raw = json.loads(pj.read_text())

        def walk(o):
            if isinstance(o, dict):
                for k, v in o.items():
                    if isinstance(v, dict):
                        walk(v)
                    elif k.endswith("_order") or k.endswith("_size"):
                        cfg[k] = int(v)

        walk(raw)
    elif pt.exists():
        for line in pt.read_text().splitlines():
            line = line.split("#")[0]
            m = re.match(r"\s*(\w+_(?:order|size))\s*=\s*(\d+)", line)
            if m:
                cfg[m.group(1)] = int(m.group(2))
    return cfg


def find_p_grid_size(app: Path, default: int = 64) -> tuple[int, str]:
    """extents[0] is integral_view.size(), i.e. the external momentum grid.

    It lives in the model header as a compile-time constant, not in the parameter file.
    """
    for name in ("model.hh", "model_vertexE.hh", "../model.hh"):
        f = app / name
        if not f.exists():
            continue
        m = re.search(r"p_grid_size\s*=\s*(\d+)", f.read_text())
        if m:
            return int(m.group(1)), f"{name}"
    return default, "DEFAULT (model.hh not found)"


def flow_integrator(flow_dir: Path) -> str | None:
    for hh in sorted(flow_dir.glob("*.hh")):
        if hh.name in ("kernel.hh", "kernels.hh", "numtrace_verdict.hh"):
            continue
        txt = hh.read_text(errors="ignore")
        for name in INTEGRATORS:
            if re.search(rf"\b{name}\s*<", txt):
                return name
    return None


def run_app(app: Path, lbs: list[int | None], nmats: int) -> None:
    cfg = read_config(app)
    pgrid, pgrid_src = find_p_grid_size(app)
    flows = sorted(p for p in (app / "flows").iterdir() if p.is_dir()) if (app / "flows").is_dir() else []

    print(f"\n=== {app}")
    print(f"    p_grid_size = {pgrid}  (from {pgrid_src});  orders: "
          + ", ".join(f"{k}={v}" for k, v in sorted(cfg.items()) if k.endswith("_order")))

    rows = []
    for fd in flows:
        integ = flow_integrator(fd)
        if integ is None:
            continue
        keys = INTEGRATORS[integ]
        missing = [k for k in keys if k != MATSUBARA and k not in cfg]
        if missing:
            print(f"    {fd.name:12s} {integ}: missing config keys {missing}")
            continue
        extents = [pgrid] + [nmats if k == MATSUBARA else cfg[k] for k in keys]
        rows.append((fd.name, integ, extents))

    if not rows:
        print("    no flows with a recognised integrator")
        return

    hdr = f"    {'flow':<12s} {'rank':>4s} {'extents':<22s}"
    for lb in lbs:
        hdr += f" {('LB=' + (str(lb) if lb else 'off')):>16s}"
    print(hdr)

    totals = {lb: [0, 0] for lb in lbs}
    for name, integ, extents in rows:
        line = f"    {name:<12s} {len(extents):>4d} {str(extents):<22s}"
        for lb in lbs:
            a = analyse(extents, lb)
            line += f" {a['block']:>4d}thr {a['waste']:>6.3f}x"
            totals[lb][0] += a["launched"]
            totals[lb][1] += a["useful"]
        print(line)

    line = f"    {'TOTAL':<12s} {'':>4s} {'':<22s}"
    for lb in lbs:
        L, U = totals[lb]
        line += f" {'':>7s} {L / U:>6.3f}x"
    print(line)

    # name the culprit dimensions at the shipped default
    print("\n    masked dimensions at the shipped LB=128:")
    any_bad = False
    for name, integ, extents in rows:
        a = analyse(extents, 128)
        if a["bad_dims"]:
            any_bad = True
            det = ", ".join(
                f"dim{i} extent {e} vs tile {t} -> {-(-e // t) * t} slots" for i, e, t in a["bad_dims"]
            )
            print(f"      {name:<12s} tile {a['tile']}  {det}")
    if not any_bad:
        print("      none - every tile dimension divides its extent")

    # what a divisibility-aware tile would give at the same 128-thread budget
    print("\n    divisibility-aware tile at a 128-thread budget (proposed):")
    print(f"      {'flow':<12s} {'kokkos tile':<20s} {'blk':>4s} {'waste':>7s}   "
          f"{'proposed tile':<20s} {'blk':>4s} {'waste':>7s}  {'gain':>7s}")
    for name, integ, extents in rows:
        a = analyse(extents, 128)
        dt = divisibility_tile(extents, 128)
        b = analyse(extents, 128, tile=dt)
        gain = a["launched"] / b["launched"]
        print(f"      {name:<12s} {str(a['tile']):<20s} {a['block']:>4d} {a['waste']:>6.3f}x   "
              f"{str(dt):<20s} {b['block']:>4d} {b['waste']:>6.3f}x  {gain:>6.3f}x")


def selftest() -> int:
    """Check the transcription against cases derivable by hand from the Kokkos source."""
    fails = 0

    def chk(label, got, want):
        nonlocal fails
        ok = got == want
        if not ok:
            fails += 1
        print(f"  [{'ok ' if ok else 'FAIL'}] {label}: got {got}, want {want}")

    # rank 5, no launch bounds: the recommended tile fits under 512 (16*2*4*2*1 = 256)
    chk("rank5 tile, LB off", kokkos_tile([64, 32, 6, 6, 6], None), [16, 2, 4, 2, 1])
    # LB=128: 16*2*4 = 128 exactly, so dim3's 2 must halve to 1
    chk("rank5 tile, LB=128", kokkos_tile([64, 32, 6, 6, 6], 128), [16, 2, 4, 1, 1])
    # LB=96: 16*2*4 = 128 > 96, so dim2 halves 4->2 giving 64; dim3 would make 128 > 96 -> 1
    chk("rank5 tile, LB=96", kokkos_tile([64, 32, 6, 6, 6], 96), [16, 2, 2, 1, 1])
    chk("rank5 tile, LB=64", kokkos_tile([64, 32, 6, 6, 6], 64), [16, 2, 2, 1, 1])
    # block size is the tile product
    chk("block at LB=128", analyse([64, 32, 6, 6, 6], 128)["block"], 128)
    chk("block at LB=96", analyse([64, 32, 6, 6, 6], 96)["block"], 64)
    # the headline numbers
    chk("nf2 waste at LB=128", round(analyse([64, 32, 6, 6, 6], 128)["waste"], 4), round(8 / 6, 4))
    chk("nf2 waste at LB=96", analyse([64, 32, 6, 6, 6], 96)["waste"], 1.0)
    chk("ym  waste at LB=128", analyse([64, 32, 8, 8, 8], 128)["waste"], 1.0)
    chk("ym  waste at LB=96", analyse([64, 32, 8, 8, 8], 96)["waste"], 1.0)
    print(f"\n{'PASS' if fails == 0 else f'{fails} FAILURES'}")
    return fails


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("apps", nargs="*", type=Path)
    ap.add_argument("--lb", default="off,64,96,128,192,256",
                    help="comma-separated launch-bounds values; 'off' means no launch bounds")
    ap.add_argument("--nmats", type=int, default=16,
                    help="assumed Matsubara axis length for finite-T flows "
                         "(runtime-determined: grid_size[dim-1] = matsubara_nodes.size())")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        return 1 if selftest() else 0
    if not args.apps:
        ap.error("give at least one app directory, or --selftest")

    lbs = [None if v.strip() == "off" else int(v) for v in args.lb.split(",")]
    for app in args.apps:
        run_app(app, lbs, args.nmats)
    return 0


if __name__ == "__main__":
    sys.exit(main())
