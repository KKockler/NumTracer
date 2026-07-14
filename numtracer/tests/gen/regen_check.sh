#!/usr/bin/env bash
# regen_check.sh — the codegen test that the rest of `ctest` is NOT.
#
# Every other test compiles the PRE-GENERATED headers in tests/gen/ and runs them; Wolfram is never
# invoked. So a green suite says nothing about a change to numtracer/mathematica/. This script
# REGENERATES the kernels from source, rebuilds, and re-runs the flow tests.
#
#   tests/gen/regen_check.sh [flow ...]   regenerate (default: the numeric flows below), rebuild, test
#   tests/gen/regen_check.sh --ctest      as above, but restore tests/gen/ on success (used by ctest)
#   tests/gen/regen_check.sh --restore    discard regenerated kernels, back to HEAD
#
# Wired into ctest as the `codegen`-labelled test `codegen_regen`:
#   cmake -B build -DNUMTRACER_TEST_CODEGEN=ON && ctest --test-dir build -L codegen
# It re-invokes ctest with `-LE codegen` so it cannot recurse into itself.
#
# Reading a failure:
#   this fails but plain ctest passes -> your codegen change altered the emitted kernel.
#   both fail                         -> the committed artifacts are stale, not necessarily wrong.
#                                        See tests/gen/README.md and FUNKIT_ROUTING_ISSUE.md.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
GEN="$REPO/numtracer/tests/gen"
BUILD="${NT_BUILD_DIR:-$REPO/numtracer/build}"

restore() { git -C "$REPO" checkout -- numtracer/tests/gen/ 2>/dev/null; }

CTEST_MODE=0
case "${1:-}" in
  --restore) restore; echo "restored tests/gen/ to HEAD"; exit 0 ;;
  --ctest)   CTEST_MODE=1; shift ;;
esac

# Every NumTracer-generated kernel a test compares. `gen_zq_collect` emits a NumTracer kernel too
# (the collected Zq path, via MakeNTKernel) even though its name does not end in `_numeric` — leaving
# it out silently left one stale kernel behind and failed flow_zq_collect.
DEFAULT_FLOWS=(
  gen_qcd_za3_numeric gen_qcd_za4_numeric gen_zaqbq1_numeric gen_zq_numeric
  gen_qcd_za_numeric gen_qcd_za3_147_numeric gen_qcd_za4_147_numeric
  gen_qcd_aqbq147_numeric gen_ym_zacbc_numeric gen_qcd_etapil_numeric
  gen_zq_collect
)
FLOWS=("$@")
[[ ${#FLOWS[@]} -eq 0 ]] && FLOWS=("${DEFAULT_FLOWS[@]}")

WOLFRAM="$(command -v wolfram || command -v wolframscript || true)"
[[ -n "$WOLFRAM" ]] || { echo "ERROR: no Wolfram kernel on PATH"; exit 1; }
[[ -d "$BUILD" ]]   || { echo "ERROR: build dir not found: $BUILD"; exit 1; }

echo "regenerating ${#FLOWS[@]} flow(s) from source ..."
fail=0
for f in "${FLOWS[@]}"; do
  [[ -f "$GEN/$f.wls" ]] || { echo "  SKIP $f (no such .wls)"; continue; }
  printf '  %-28s ' "$f"
  if (cd "$GEN" && "$WOLFRAM" -script "$f.wls" >"/tmp/regen_$f.log" 2>&1); then
    echo "ok"
  else
    echo "FAILED (see /tmp/regen_$f.log)"; fail=1
  fi
done
if (( fail )); then
  echo "generation failed; leaving tests/gen/ as-is for inspection (restore: $0 --restore)"
  exit 1
fi

echo "rebuilding against the regenerated kernels ..."
if ! cmake --build "$BUILD" -j "$(nproc)" >/tmp/regen_build.log 2>&1; then
  echo "ERROR: build failed (see /tmp/regen_build.log)"; exit 1
fi

echo "running the flow tests against the REGENERATED kernels ..."
# -LE codegen: never re-enter this test (it is itself a ctest when NUMTRACER_TEST_CODEGEN=ON)
if ctest --test-dir "$BUILD" -LE codegen --output-on-failure; then
  echo
  echo "PASS — codegen reproduces kernels that still match the oracles."
  if (( CTEST_MODE )); then
    restore
    # the test binaries now embed the regenerated kernels; rebuild so the tree is left consistent
    cmake --build "$BUILD" -j "$(nproc)" >/dev/null 2>&1
    echo "(restored tests/gen/ to HEAD)"
  fi
  exit 0
else
  echo
  echo "FAIL — the regenerated kernels differ from the frozen oracles."
  echo "Before blaming Codegen.m, check whether the ORACLE is the stale side:"
  echo "  * the refshim FORM oracles are frozen MANUAL copies of flows/<NAME>/kernel.hh (README.md)"
  echo "  * a FunKit backend/routing change moves the integrand without changing the physics"
  echo "    (FUNKIT_ROUTING_ISSUE.md at the repo root)"
  echo "tests/gen/ left as-is for inspection; restore with: $0 --restore"
  exit 1
fi
