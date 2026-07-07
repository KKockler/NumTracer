/// @file export.hpp
/// @brief Header-only vs. separately-compiled toggle for the build-time engine/generator.
///
/// NumTracer ships, by default, as a compiled static library (`NumTracer::NumTracer`): the heavy
/// engine bodies (numeric contraction, SU(N) folding, lowering/CSE, polynomial arithmetic) are
/// declared in the headers but *defined* once in `src/numtracer.cpp`, so a per-flow generator TU —
/// or any consumer — only parses the declarations and links the prebuilt object instead of
/// re-instantiating and `-O2`-optimising the whole engine on every translation unit.
///
/// Define `NUMTRACER_HEADER_ONLY=1` (the `NumTracer::NumTracer_headeronly` target) to fall back to
/// the classic header-only build, where every engine function is `inline` and self-contained — for
/// pure-header or cross-toolchain use where linking a local `.a` is undesirable.
///
/// This toggle affects only how the *engine/generator* is compiled. The kernel the generator emits
/// is unaffected: it is self-contained straight-line C++ against `codegen/runtime.hpp` and never
/// includes these headers.
///
/// A split engine header keeps forward declarations of its public entry points always visible, marks
/// their definitions @ref NUMTRACER_FUNC, and wraps its whole definition region (public entry points +
/// internal helpers) in `#if NUMTRACER_DEFINE_BODIES`. That region compiles in exactly two situations:
/// in a header-only build, or inside the single library TU `src/numtracer.cpp` (which defines
/// `NUMTRACER_COMPILING_LIB`). In a normal compiled consumer TU the region is skipped — only the
/// declarations are seen, and the definitions are linked from `NumTracer::NumTracer`.
#pragma once

#ifndef NUMTRACER_HEADER_ONLY
#define NUMTRACER_HEADER_ONLY 0
#endif

/// @brief True in the TUs that must emit the engine bodies: header-only builds, or the library TU.
#if NUMTRACER_HEADER_ONLY || defined(NUMTRACER_COMPILING_LIB)
#define NUMTRACER_DEFINE_BODIES 1
#else
#define NUMTRACER_DEFINE_BODIES 0
#endif

#if NUMTRACER_HEADER_ONLY
/// @brief Definition-site marker for a split engine function: `inline` in header-only mode.
#define NUMTRACER_FUNC inline
#else
/// @brief Definition-site marker for a split engine function: external linkage (compiled once).
#define NUMTRACER_FUNC
#endif
