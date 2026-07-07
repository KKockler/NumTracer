/// @file core/config.hpp
/// @brief Centralised compile-time build tunables and portable compiler-pragma shims.
///
/// Gathered here so the trade-offs live in one place rather than being scattered through the
/// headers that use them. Nothing here changes *results* — only the build's cost profile or a
/// compiler-specific workaround.
#pragma once

#include <cstddef>

// ---- portable compiler pragmas -------------------------------------------------------------
// Region guard that disables loop auto-vectorization for the enclosed function definitions.
// Sole use: working around a GCC 16 -O3 -march=native (AVX-512VL) miscompile of the complex-`Cx`
// std::vector fold loops in dense/dtensor.hpp (see there). Clang -O3 -march=native is NOT affected,
// so its branch is a documented no-op; any other compiler gets the default no-op. Bracket the
// affected definitions with `NT_BEGIN_NO_LOOP_VECTORIZE` ... `NT_END_NO_LOOP_VECTORIZE`.
#if defined(__GNUC__) && !defined(__clang__)
#define NT_BEGIN_NO_LOOP_VECTORIZE _Pragma("GCC push_options") _Pragma("GCC optimize(\"no-tree-loop-vectorize\")")
#define NT_END_NO_LOOP_VECTORIZE _Pragma("GCC pop_options")
#elif defined(__clang__)
#define NT_BEGIN_NO_LOOP_VECTORIZE // clang does not exhibit the GCC AVX-512VL fold bug
#define NT_END_NO_LOOP_VECTORIZE
#else
#define NT_BEGIN_NO_LOOP_VECTORIZE // unknown compiler: no-op
#define NT_END_NO_LOOP_VECTORIZE
#endif