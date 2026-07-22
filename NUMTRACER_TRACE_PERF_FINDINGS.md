# Large-trace generator efficiency — findings

Investigation of why the dense 4-point quark-gluon flows (ZAAqbq1/ZAAqbq2) cost so much generator
RAM (~25 GB unthrottled) and time. Diagnosed 2026-07-22 across five angles: single-contraction RAM,
per-contraction runtime (phase A), fold runtime (phase B), Wolfram-side trace-set size, and the
polynomial data structure / allocator. Every number below was measured (micro-benchmarks against the
real headers, or `sizeof` probes); the per-angle detail and unified-diff sketches are archived with
the investigation.

Baselines (from `NUMTRACER_LOWERING_RAM.md` / `NUMTRACER_PERF_NEXT.md`):
- ZAAqbq1 generator ~53 s (≈27 s phase A + 4 s phase B + compile), single contraction >8 GB peak.
- ZAAqbq2 generator ~206 s (≈165 s phase A + 40.5 s phase B).
- Polynomial term (`pair<Mono,Cx>`) = **128 B measured**; `poly_bytes = t.size()·128`. So the 41 MB
  net poly ≈ 320k monomials and the >8 GB single contraction ≈ 63M live terms.

## Two insights that connect the angles

1. **One `Mono` change serves the top RAM finding *and* the top per-contraction runtime finding.**
   The RAM angle wants exponents bit-packed to shrink the term (128→72 B, −44%; >8 GB → ~4.5 GB). The
   phase-A angle wants a packed `uint64` exponent *key* to replace the element-wise `Mono::operator<`
   walk in `from_scratch`'s sort (measured −7…9%). These are the same packed-exponent representation.
2. **A `tcmalloc`/`jemalloc` preload is a zero-code lever both runtime angles surfaced
   independently** (−7% from tiny-alloc churn; also lowers the fragmentation floor). One build flag on
   the disposable generator binary.

## Root cause of the size

The `{A,A,qb,q}` flow's internal quark-gluon vertices each carry 3 tensor structures, so a diagram
with n such vertices distributes as ~3ⁿ (12,132 contractions for ZAAqbq1). That is physics — but the
**distribution is not forced**: the codebase already has an eager "dressing collection" path that
keeps a structure sum as one contracted trace (used for propagator numerators `Mq·δ − iZ·p̸`), and it
is only *gated off* for the multi-structure vertex because the eager form cannot yet carry the
vertex's free gluon leg (`collectibleDiracSumQ`, `DSL.m:221-234`; C++ `DSlotOpt`,
`numeric_contract.hpp:943-951`). Enabling it is the exact, physics-preserving version of the
`dropAqbq47` approximation.

## Ranked levers

| # | Lever | Axis | Measured / estimated impact | Risk |
|---|---|---|---|---|
| **B1** | `MPolyFactory::scaled` — kill the dead `std::sort` in `scale_trace`/`scaleCx` (constant scaling currently routed through full `operator*`) | phase-B runtime | 2.8–6.3× on the op → **~15–30% off phase B** (ZAAqbq2 40.5→~28–34 s) | near-zero — bit-identical, per-leaf, no reassociation |
| **T** | `tcmalloc_minimal` preload on the generator binary | runtime + RAM floor | −7% runtime; less fragmentation | zero — disposable binary |
| **D2** | 32-bit `size_type` on the two `Mono` small_vectors | RAM | term 128→112 (−12.5%); may clear the 10 GB cap alone | near-zero — 2 `using` lines, no arithmetic change |
| **A1+D1** | Packed-exponent `Mono` (uint64 sort key + packed storage, exact `e`-fallback for degree>31 / nsym>12) | RAM **and** runtime | term 128→72 (−44%), >8 GB → ~4.5 GB, **and** −7…9% per contraction | medium — touches every `e[k]` site; both agents kept an exact fallback |
| **R1** | Threshold-gated **blocked `operator*`** — cap transient scratch above ~2²⁰; hot path (mean-4-term) byte-identical | single-contraction RAM | collapsing >8 GB spike → **~0.3–0.5 GB (15–25×)** | low — block-sum reassociation < the shipped phase-B tree fold |
| **W1** | Route the internal quark-gluon vertex through **eager dressing collection** instead of 3ⁿ distribution | generator RAM (root cause) + **retires `dropAqbq47`** | ~15–20× fewer emitted nets (488→~24), table 20.1→~1–2 MB | low for emitted-size (exact); the contraction-RAM part needs a matrix-of-`DPoly` follow-on |

**Composition.** A1+D1 (fewer bytes/term) × R1 (fewer live transient terms) multiply on the single
contraction. W1 is the only lever that removes the 3ⁿ root cause *exactly* — the physics payoff (full
ZAqbq1/4/7 basis, no `dropAqbq47`); its emitted-size part is self-contained, the contraction-RAM part
wants the transfer-matrix (matrix-of-`DPoly`) contraction as a follow-on.

## Measured and rejected — do not chase

- **Constant×poly fast-path** skipping the sort in `operator*`: measured ~0 (const-operand products
  come out already sorted; insertion-sort is O(n)-free there — confirms the "at its floor" note).
- **Narrowing the term 128→96 B** for sort speed: measured ~0 — the sort is **compare-bound, not
  move-bound**, which is exactly what points to the packed-key lever A1.
- **Phase-B in-place `operator+=`**: only ~5–7% (the merge copies are inherent; O(log n) partials live
  at once forbid a single ping-pong buffer).
- **Real (non-complex) `Cx`**: ~6% RAM only (term 128→120). Worth doing for compute / vacuum-Re
  reasons, not as a RAM lever.
- Everything already on the pre-existing rejected list (`MPoly::t` inline, matmul sparsity skip,
  hash-collect in `from_scratch`, single-pass `Mono::operator<`).

See `NUMTRACER_TRACE_PERF_PLAN.md` for the staged, measured implementation plan.
