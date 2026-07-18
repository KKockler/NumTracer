cmake_minimum_required(VERSION 3.20)
project(NumTracer VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Default to a per-user prefix (same convention as DiFfRG: consumers say
# `find_package(NumTracer REQUIRED HINTS ~/.local/share/NumTracer)`), but never
# override an explicit -DCMAKE_INSTALL_PREFIX.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/.local/share/NumTracer"
      CACHE PATH "NumTracer install prefix" FORCE)
endif()

# Deliberately NO fast-math flags. `-fno-math-errno` measured ~8% faster on the A4 kernel
# (752 vs 814 ns) and is itself NaN-safe (only skips errno on libm calls), but `-ffast-math`
# / `-ffinite-math-only` assume no NaN/Inf and silently mask them — unacceptable for fRG
# kernels where a NaN signals a real problem. Keep the build NaN-honest; the win is small and
# also speeds the FORM reference, so it does not close the NumTracer-vs-FORM gap anyway.
# CXX-only: once the GPU tests enable_language(CUDA), these must not leak into nvcc
# (-march=native is not an nvcc flag; CUDA targets get their own flags below).
# Directory-scoped, so nothing of it is exported to consumers of the installed package.
add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-O3;-march=native;-Wall;-Wextra>")

# Opt-in sanitizers for validation builds (off by default; never affects the shipped
# package). Accepts a semicolon list of sanitizer names, e.g.
#   cmake -S . -B build-asan -DNUMTRACER_SANITIZE="ADDRESS;UNDEFINED"
# Recognised: ADDRESS, UNDEFINED, THREAD, MEMORY (case-insensitive). When enabled we drop
# to -O1 and keep frame pointers so the reports carry clean, inlining-free stack traces;
# the flags go to both compile and link lines (sanitizer runtimes must be linked in).
set(NUMTRACER_SANITIZE "OFF" CACHE STRING
    "Sanitizers for CXX builds: OFF, or a list of ADDRESS;UNDEFINED;THREAD;MEMORY")
if(NOT NUMTRACER_SANITIZE STREQUAL "OFF")
  set(_nt_san_list "")
  foreach(_s IN LISTS NUMTRACER_SANITIZE)
    string(TOLOWER "${_s}" _s_lc)
    if(_s_lc STREQUAL "address")
      list(APPEND _nt_san_list "address")
    elseif(_s_lc STREQUAL "undefined")
      list(APPEND _nt_san_list "undefined")
    elseif(_s_lc STREQUAL "thread")
      list(APPEND _nt_san_list "thread")
    elseif(_s_lc STREQUAL "memory")
      list(APPEND _nt_san_list "memory")
    else()
      message(FATAL_ERROR "NUMTRACER_SANITIZE: unknown sanitizer '${_s}'")
    endif()
  endforeach()
  string(REPLACE ";" "," _nt_san "${_nt_san_list}")
  message(STATUS "NumTracer: sanitizers enabled (-fsanitize=${_nt_san}, -O1)")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-O1;-g;-fno-omit-frame-pointer;-fsanitize=${_nt_san}>")
  add_link_options("-fsanitize=${_nt_san}")
endif()

# ##############################################################################
# The library: a compiled static library by default, header-only on opt-in
# ##############################################################################
#
# The build-time engine/generator (numeric contraction, SU(N) folding, lowering/CSE, polynomial
# arithmetic) is heavy, template/constexpr-dense header code. Compiled once into this static library
# its bodies no longer re-instantiate and -O2-optimise on every consumer/generator TU — which is what
# made every generation run pay a multi-second compile floor regardless of kernel size. The engine
# functions are declared in the headers and defined once in src/numtracer.cpp (see core/export.hpp).
# This affects only how the engine ships; the kernel it *emits* stays standalone header-only C++.
add_library(NumTracer STATIC src/numtracer.cpp)
add_library(NumTracer::NumTracer ALIAS NumTracer)
target_compile_features(NumTracer PUBLIC cxx_std_20)
target_include_directories(NumTracer PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

# Header-only fallback: same headers, but every engine function is `inline` and self-contained, so no
# `.a` link is needed. For pure-header or cross-toolchain use. A given binary must pick ONE mode across
# all its TUs (ODR); linking NumTracer::NumTracer is the default, this is the escape hatch.
add_library(numtracer_headeronly INTERFACE)
set_target_properties(numtracer_headeronly PROPERTIES EXPORT_NAME NumTracer_headeronly)
add_library(NumTracer::NumTracer_headeronly ALIAS numtracer_headeronly)
target_compile_features(numtracer_headeronly INTERFACE cxx_std_20)
target_compile_definitions(numtracer_headeronly INTERFACE NUMTRACER_HEADER_ONLY=1)
target_include_directories(numtracer_headeronly INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

# ##############################################################################
# Install & export (headers + CMake package config; arch-independent)
# ##############################################################################

include(CMakePackageConfigHelpers)

install(DIRECTORY include/ DESTINATION include)
install(TARGETS NumTracer numtracer_headeronly EXPORT NumTracerTargets)
install(EXPORT NumTracerTargets NAMESPACE NumTracer::
        DESTINATION lib/cmake/NumTracer)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/NumTracerConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/NumTracerConfig.cmake
  INSTALL_DESTINATION lib/cmake/NumTracer)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/NumTracerConfigVersion.cmake
  COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/NumTracerConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/NumTracerConfigVersion.cmake
        DESTINATION lib/cmake/NumTracer)

# ##############################################################################
# Mathematica package install (to $UserBaseDirectory/Applications, FunKit-style)
# ##############################################################################
#
# The Wolfram front-end (NumTrace / MakeNTKernel) installs beside FunKit so that
# `Needs["NumTracer`"]` works from anywhere. The Applications dir is queried from
# wolframscript at configure time; if no Wolfram installation is found the install
# is skipped with a warning (the C++ library install is unaffected).
option(NUMTRACER_INSTALL_MATHEMATICA
       "Install the Wolfram package to $UserBaseDirectory/Applications" ON)

if(NUMTRACER_INSTALL_MATHEMATICA)
  find_program(WolframLanguage_wolframscript_EXE NAMES wolframscript wolfram math)
  if(DEFINED ENV{WOLFRAM_SCRIPT})
    set(WolframLanguage_wolframscript_EXE $ENV{WOLFRAM_SCRIPT})
  endif()
  if(WolframLanguage_wolframscript_EXE)
    execute_process(
      COMMAND ${WolframLanguage_wolframscript_EXE} -script
              ${CMAKE_CURRENT_SOURCE_DIR}/cmake/get_wolfram_app_dir.m
      OUTPUT_VARIABLE WOLFRAM_APP_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE _wolfram_rc)
    string(REPLACE "\n" "" WOLFRAM_APP_DIR "${WOLFRAM_APP_DIR}")
    if(_wolfram_rc EQUAL 0 AND WOLFRAM_APP_DIR)
      message(STATUS "Wolfram Language application directory: ${WOLFRAM_APP_DIR}")
      # Record the header install location for the package: the "Generate" backend
      # compiles the emitted generator with g++ against these headers.
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NumTracerPaths.m.in
                     ${CMAKE_CURRENT_BINARY_DIR}/NumTracerPaths.m @ONLY)
      install(FILES mathematica/NumTracer.m
                    mathematica/DSL.m
                    mathematica/Frames.m
                    mathematica/Codegen.m
                    mathematica/FunKitAdapter.m
                    mathematica/DiFfRG_compat.m
                    mathematica/PacletInfo.m
                    ${CMAKE_CURRENT_BINARY_DIR}/NumTracerPaths.m
              DESTINATION ${WOLFRAM_APP_DIR}/NumTracer)
    else()
      message(WARNING "wolframscript found but querying $UserBaseDirectory/Applications "
                      "failed — skipping the Mathematica package install")
    endif()
  else()
    message(WARNING "wolframscript not found — skipping the Mathematica package install "
                    "(set the WOLFRAM_SCRIPT environment variable, or configure with "
                    "-DNUMTRACER_INSTALL_MATHEMATICA=OFF to silence this warning)")
  endif()
endif()

# ##############################################################################
# Tests & benchmarks (top-level builds only; off under add_subdirectory)
# ##############################################################################

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(NT_IS_TOP_LEVEL ON)
else()
  set(NT_IS_TOP_LEVEL OFF)
endif()
option(NUMTRACER_BUILD_TESTS "Build the NumTracer tests and benchmarks" ${NT_IS_TOP_LEVEL})

if(NUMTRACER_BUILD_TESTS)
  # The tests compile committed machine-generated kernels (tests/gen/, tests/refshim/). NumTracer's own
  # emitter now marks the unused-by-construction bits [[maybe_unused]] (see codegen/gen.hpp and
  # mathematica/Codegen.m), so freshly regenerated NumTracer kernels are -Wunused clean. This fallback
  # covers the two remaining sources that are NOT fixable at build time: committed NumTracer kernels not
  # yet regenerated against the current emitter, and the FormTracer/FunKit-emitted oracle kernels
  # (external — we don't own their emitter). Silence just these two warnings for the internal test/bench
  # targets only — directory-scoped and CXX-only, so it never reaches the shipped library (added above),
  # the installed package, or nvcc; -Wall -Wextra otherwise stays in force.
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-Wno-unused-parameter;-Wno-unused-variable>")

  # On-the-fly color-contraction.
  add_executable(sun_contract tests/sun_contract.cpp)
  target_link_libraries(sun_contract PRIVATE NumTracer::NumTracer)

  # Numeric (matrix-product) contraction backend (numeric/*): MPoly-entry gamma trace via 4x4
  # matrix products + projector index-sum, vs chiral_gamma_trace and direct numeric contraction;
  # validates the monomial-vs-non-monomial projector-atom cancellation. The Part-B correctness gate.
  add_executable(test_numeric_contract tests/test_numeric_contract.cpp)
  target_link_libraries(test_numeric_contract PRIVATE NumTracer::NumTracer)

  # Projector identities for the numeric backend: transverse P_T and longitudinal P_L factors
  # (nprojT / nprojL) checked against the 4x4 oracle and the algebraic identities (P_T+P_L=δ,
  # idempotency, P_T·P_L=0, traces 3/1, transversality/longitudinality).
  add_executable(test_projector_identity tests/test_projector_identity.cpp)
  target_link_libraries(test_projector_identity PRIVATE NumTracer::NumTracer)

  # Projector-algebra fusion (ndetail::fuse_projectors): same-momentum projectors sharing a dummy
  # index collapse pre-expansion — idempotency P·P→P / tr P, orthogonality L⟂{T,E,M}, E⟂M → 0;
  # structural folds (all symmetric slot variants, chains, traces, no-op negatives) + value-preserving.
  add_executable(test_projector_fusion tests/test_projector_fusion.cpp)
  target_link_libraries(test_projector_fusion PRIVATE NumTracer::NumTracer)

  # Dressing-collection layer (numeric/dpoly.hpp + numeric_value_dressed): DPoly arithmetic, the
  # dressing-free reduction (DPoly == MPoly path, same emitted program), and a dressed quark line
  # collected-vs-distributed. The symbolic-dressing-collection correctness gate.
  add_executable(test_dpoly tests/test_dpoly.cpp)
  target_link_libraries(test_dpoly PRIVATE NumTracer::NumTracer)

  # The build-time numeric colour/flavour net (network/sun_net.hpp): fundamental T^a, two groups,
  # vs analytic group factors.
  add_executable(test_col_numeric tests/test_col_numeric.cpp)
  target_link_libraries(test_col_numeric PRIVATE NumTracer::NumTracer)

  # THE SUB-TERM DEDUP GATE — the only flow with real sub-term density (~7.5k contractions over ~820
  # distinct traces, 9x redundancy, longest per-net fold ~140). Every other committed flow has exactly
  # ONE sub-term per net, so the dedup's merge / trace cache / tree fold are no-ops there. Grades the
  # dedup-ON kernel against a dedup-OFF (NT_GEN_NO_DEDUP) control generated from the same derivation.
  add_executable(compare_zaaqbq1_small tests/refshim/compare_zaaqbq1_small.cpp)
  target_include_directories(compare_zaaqbq1_small PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_zaaqbq1_small PRIVATE NumTracer::NumTracer)

  # Silent-corruption guards on paths no committed flow exercises: network::add refusing a complex
  # vector-sum coefficient (real-only Elem::vlc), and the empty colour net folding to the identity.
  add_executable(test_guards tests/test_guards.cpp)
  target_link_libraries(test_guards PRIVATE NumTracer::NumTracer)

  # The generator's sub-term dedup (numeric/trace_fold.hpp): the balanced tree fold vs the left fold
  # it replaces, and every trace-cache regime. This is the ONLY test that sees the fold: the dedup
  # needs many sub-terms per net, and every flow in tests/gen has exactly one — so the flow tests and
  # even the codegen test are blind to it. Includes a work bound (a botched carry still returns the
  # right VALUE while quietly degenerating back into the left fold, which no value check can see).
  add_executable(test_trace_fold tests/test_trace_fold.cpp)
  target_link_libraries(test_trace_fold PRIVATE NumTracer::NumTracer)

  # Phase-2 example 2: pion <- quark loop -> pion (Yukawa) — fundamental flavour x Dirac x colour.
  add_executable(compare_pion_quark tests/refshim/compare_pion_quark.cpp)
  target_include_directories(compare_pion_quark PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_pion_quark PRIVATE NumTracer::NumTracer)

  # Phase-2 example 3: gluon <- quark loop -> gluon — fundamental colour via the numeric path.
  add_executable(compare_glu_quark tests/refshim/compare_glu_quark.cpp)
  target_include_directories(compare_glu_quark PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_glu_quark PRIVATE NumTracer::NumTracer)

  add_executable(compare_zaqbq1 tests/refshim/compare_zaqbq1.cpp)
  target_include_directories(compare_zaqbq1 PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_zaqbq1 PRIVATE NumTracer::NumTracer)

  # Symbolic dressing collection: collected Zq (DressingCollection -> True) vs the distributed Zq.
  add_executable(compare_zq_collect tests/refshim/compare_zq_collect.cpp)
  target_include_directories(compare_zq_collect PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_zq_collect PRIVATE NumTracer::NumTracer)

  # NUMERIC (matrix-product) backend validation: ZA gluon self-energy vs the FORM oracle; ZA3 quark
  # triangle vs the FORM oracle (the headline — generates in seconds at <0.4 GB where the inv path OOMs).
  add_executable(compare_za_num tests/refshim/compare_za_num.cpp)
  target_include_directories(compare_za_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_za_num PRIVATE NumTracer::NumTracer)

  add_executable(compare_za3_num tests/refshim/compare_za3_num.cpp)
  target_include_directories(compare_za3_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_za3_num PRIVATE NumTracer::NumTracer)

  # pure-YM ghost-gluon vertex: numeric vs the FORM oracle. Guards the disconnected-Lorentz-component
  # (product-not-sum) codegen bug — the only flow with such a diagram.
  add_executable(compare_zacbc_num tests/refshim/compare_zacbc_num.cpp)
  target_include_directories(compare_zacbc_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_zacbc_num PRIVATE NumTracer::NumTracer)

  # ZA3 with the FULL quark-gluon vertex {1,4,7} tensor basis (the dense struct-7 case). Gate for the
  # struct-7 codegen fix: with ZAqbq4=ZAqbq7=0 the 1/4/7 numeric kernel must reduce EXACTLY to the
  # FORM-validated struct-1 numeric kernel (ZA3_num_kernel), and differ from it with the real 4/7
  # dressings. (A full oracle is impractical at this size; the struct-1 kernel is the trusted reference.)
  add_executable(compare_za3_147_num tests/refshim/compare_za3_147_num.cpp)
  target_include_directories(compare_za3_147_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_za3_147_num PRIVATE NumTracer::NumTracer)

  # ZA3 1/4/7 numeric vs the copied FormTracer oracle (independent physical check of structures 4/7).
  add_executable(compare_za3_147_form tests/refshim/compare_za3_147_form.cpp)
  target_include_directories(compare_za3_147_form PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_za3_147_form PRIVATE NumTracer::NumTracer)

  # etaPiL pion anomalous dimension: numeric kernel vs the copied FormTracer oracle (independent
  # tensor-trace engines on the same pion two-point wave-function loop).
  add_executable(compare_etapil_form tests/refshim/compare_etapil_form.cpp)
  target_include_directories(compare_etapil_form PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_etapil_form PRIVATE NumTracer::NumTracer)

  # ns/eval timing of the 1/4/7 numeric kernel vs the FORM oracle (build on demand: not a pass/fail
  # test). `cmake --build build --target bench_147 && ./build/bench_147` reproduces the README numbers.
  add_executable(bench_147 EXCLUDE_FROM_ALL tests/refshim/bench_147.cpp)
  target_include_directories(bench_147 PRIVATE tests/gen tests/refshim)
  target_link_libraries(bench_147 PRIVATE NumTracer::NumTracer)

  add_executable(compare_za4_num tests/refshim/compare_za4_num.cpp)
  target_include_directories(compare_za4_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_za4_num PRIVATE NumTracer::NumTracer)

  # Pure-GAUGE ZA4: the topologies compare_za4_num.cpp does NOT cover (it is the QCD quark box).
  # Grades the ZA3^4 box, ZA3^2.ZA4 triangle, ZA4^2 bubble and ghost box SEPARATELY by polynomial
  # subtraction — a full-strength compare is dominated by the box and hides the others.
  add_executable(compare_ym_za4_num tests/refshim/compare_ym_za4_num.cpp)
  target_include_directories(compare_ym_za4_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_ym_za4_num PRIVATE NumTracer::NumTracer)

  # ZA4 (four-gluon box, quark loop) with the FULL {1,4,7} vertex — the heaviest 1/4/7 target, made
  # tractable by the σ^{μν} commutator fold (~14 min generation vs the multi-hour unfolded contraction).
  # The committed kernel is ~5.3 MB (tests/gen/ZA4_147_num_kernel{,s}.hh) and compiles heavily at -O3
  # (~13 s, ~0.4 GB), so this test is OFF by default rather than absent. Regenerate the kernel with
  #   wolfram -script tests/gen/gen_qcd_za4_147_numeric.wls
  # then configure -DNUMTRACER_BUILD_ZA4_147=ON (the EXISTS guard below also covers a pruned kernel).
  # Validation is by struct-reduction: with ZAqbq4=ZAqbq7=0 the 1/4/7 kernel must reduce to the
  # FORM-validated ZA4_num_kernel (measured max rel err 8.3e-12), and differ with the real 4/7 dressings.
  option(NUMTRACER_BUILD_ZA4_147 "Build the ZA4 1/4/7 struct-reduction test (heavy ~5.3MB committed kernel)" OFF)
  if(NUMTRACER_BUILD_ZA4_147 AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/gen/ZA4_147_num_kernel.hh")
    add_executable(compare_za4_147_num tests/refshim/compare_za4_147_num.cpp)
    target_include_directories(compare_za4_147_num PRIVATE tests/gen tests/refshim)
    target_link_libraries(compare_za4_147_num PRIVATE NumTracer::NumTracer)
    add_test(NAME flow_za4_147_num COMMAND compare_za4_147_num)
  endif()

  add_executable(compare_zq_num tests/refshim/compare_zq_num.cpp)
  target_include_directories(compare_zq_num PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_zq_num PRIVATE NumTracer::NumTracer)

  # Per-flavour dressings: dressing the u-quark and d-quark DIFFERENTLY (broken flavour symmetry).
  # Two explicit flavour loops, each with its own runtime dressing (Gu/Gd, GuDot/GdDot) and NO
  # SU(Nf) flavour group; numeric kernel vs the FORM oracle + the per-flavour independence/breaking
  # checks (kernels generated by tests/gen/gen_flavour_split.wls).
  add_executable(compare_flavour_split tests/refshim/compare_flavour_split.cpp)
  target_include_directories(compare_flavour_split PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_flavour_split PRIVATE NumTracer::NumTracer)

  # Per-component (group-diagonal) dressings: dress u/d differently WITHIN the SU(2) flavour group
  # (no per-flavour diagram split). Validated by exact equivalence to flow_flavour_split (the
  # separate-diagram kernel), so no external oracle is needed (kernels by gen/gen_flavour_ingroup.wls).
  add_executable(compare_flavour_ingroup tests/refshim/compare_flavour_ingroup.cpp)
  target_include_directories(compare_flavour_ingroup PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_flavour_ingroup PRIVATE NumTracer::NumTracer)

  # Per-component ADJOINT dressing (gluon condensation): an SU(3)-adjoint line dressed component-by-
  # component (Z_A^a, 8 components) folded through the group. Validated by the collapse property
  # (all components equal == the colour-blind kernel) + per-component breaking — no external oracle
  # needed (kernels by gen/gen_gluon_condensate.wls).
  add_executable(compare_gluon_condensate tests/refshim/compare_gluon_condensate.cpp)
  target_include_directories(compare_gluon_condensate PRIVATE tests/gen tests/refshim)
  target_link_libraries(compare_gluon_condensate PRIVATE NumTracer::NumTracer)

  # Longitudinal-projector codegen gate: numeric backend vs a direct analytic oracle for a hand-built
  # net mixing ntLongProj + ntTransProj (kernel generated by tests/gen/gen_proj_numeric.wls).
  add_executable(compare_proj_num tests/compare_proj_num.cpp)
  target_include_directories(compare_proj_num PRIVATE tests/gen)
  target_link_libraries(compare_proj_num PRIVATE NumTracer::NumTracer)

  # Multi-term (Plus) Dirac projector codegen gate: numeric backend vs a direct analytic oracle for a
  # hand-built net whose projector is an EAGER Plus of MIXED gamma parity (branches with 0 and 1 γ),
  # kernel generated by tests/gen/gen_multiterm_proj_numeric.wls. The regression guard for the
  # odd-trace filter in DSL.m: judging parity on the DIAGRAM-GLOBAL γ count discarded the whole
  # diagram (global count 3 → odd) along with its non-vanishing even branch, silently emitting an
  # identically-zero kernel. The verdict must be per-branch — drop only when EVERY branch is odd.
  add_executable(compare_multiterm_proj_num tests/compare_multiterm_proj_num.cpp)
  target_include_directories(compare_multiterm_proj_num PRIVATE tests/gen)
  target_link_libraries(compare_multiterm_proj_num PRIVATE NumTracer::NumTracer)

  # Finite-temperature codegen gate: numeric backend vs a direct analytic oracle for a hand-built
  # finite-T net exercising ntElectricProj/ntMagneticProj/ntTransProj + ntSPS + ntVec[q,0] over
  # propFrameFT (kernel generated by tests/gen/gen_ftproj_numeric.wls).
  add_executable(compare_ftproj_num tests/compare_ftproj_num.cpp)
  target_include_directories(compare_ftproj_num PRIVATE tests/gen)
  target_link_libraries(compare_ftproj_num PRIVATE NumTracer::NumTracer)

  # Disconnected Dirac x Lorentz codegen gate: numeric backend vs a direct analytic oracle for a
  # hand-built net whose diagrams factorize into a closed Dirac trace, a constant colour trace, and an
  # index-disjoint pure-Lorentz scalar that MULTIPLY (kernel generated by tests/gen/gen_discdirac_numeric.wls).
  # The regression guard for MakeNTKernel::disconnectmix (the hSigL sigma-quark vertex): the numeric
  # backend factors the disconnected Lorentz scalar as its own trace and multiplies it in.
  add_executable(compare_discdirac_num tests/compare_discdirac_num.cpp)
  target_include_directories(compare_discdirac_num PRIVATE tests/gen)
  target_link_libraries(compare_discdirac_num PRIVATE NumTracer::NumTracer)

  # σL δ-closed meson-loop regression: a quark line threading two pion γ5 vertices closed by a BARE
  # spinor δ external projector (the σL scalar component) — tr(δ·γ5·S·γ5·S) ∝ l1·p1. Guards the
  # orderDiracFacs undirected-walk fix: the old directed walk collapsed this δ-closed loop to a single
  # γ5 → tr(γ5)=0. Numeric (gamma-trace theorem) vs a direct analytic oracle for the closed-form trace
  # (l1·p1)·G(l1²)·G(p²), which is non-vanishing (kernel generated by tests/gen/gen_sigl_dirac.wls).
  add_executable(compare_sigl_dirac tests/refshim/compare_sigl_dirac.cpp)
  target_include_directories(compare_sigl_dirac PRIVATE tests/gen)
  target_link_libraries(compare_sigl_dirac PRIVATE NumTracer::NumTracer)

  # Finite-temperature quark self-energy (hand-coded, dependency-free): the numeric engine contracts
  # the closed Dirac trace tr[(γ·p) γ^μ (γ·q) γ^ν] against the finite-T gluon line
  # ZAE·P_E(l) + ZAM·P_M(l) in the heat-bath frame (p0 = πT), checked against the closed-form trace
  # identity oracle. (The codegen DENSE oracle overflows the stack at -O3 for a 4-γ × projector trace;
  # this validates the same physics through the numeric path without that fragility.)
  add_executable(test_finiteT_quark tests/test_finiteT_quark.cpp)
  target_link_libraries(test_finiteT_quark PRIVATE NumTracer::NumTracer)

  enable_testing()
  add_test(NAME numeric_contract COMMAND test_numeric_contract)
  add_test(NAME dpoly COMMAND test_dpoly)
  add_test(NAME flow_zq_collect COMMAND compare_zq_collect)
  add_test(NAME projector_identity COMMAND test_projector_identity)
  add_test(NAME projector_fusion COMMAND test_projector_fusion)
  add_test(NAME proj_num COMMAND compare_proj_num)
  add_test(NAME multiterm_proj_num COMMAND compare_multiterm_proj_num)
  add_test(NAME ftproj_num COMMAND compare_ftproj_num)
  add_test(NAME discdirac_num COMMAND compare_discdirac_num)
  add_test(NAME flow_sigl_dirac COMMAND compare_sigl_dirac)
  add_test(NAME finiteT_quark COMMAND test_finiteT_quark)
  add_test(NAME col_numeric COMMAND test_col_numeric)
  add_test(NAME guards COMMAND test_guards)
  add_test(NAME trace_fold COMMAND test_trace_fold)
  add_test(NAME flow_zaaqbq1_small COMMAND compare_zaaqbq1_small)
  add_test(NAME flow_zaqbq1 COMMAND compare_zaqbq1)
  add_test(NAME flow_pion_quark COMMAND compare_pion_quark)
  add_test(NAME flow_glu_quark COMMAND compare_glu_quark)
  add_test(NAME flow_za_num COMMAND compare_za_num)
  add_test(NAME flow_za3_num COMMAND compare_za3_num)
  add_test(NAME flow_zacbc_num COMMAND compare_zacbc_num)
  add_test(NAME flow_za3_147_num COMMAND compare_za3_147_num)
  add_test(NAME flow_za3_147_form COMMAND compare_za3_147_form)
  add_test(NAME flow_etapil_form COMMAND compare_etapil_form)
  add_test(NAME flow_za4_num COMMAND compare_za4_num)
  add_test(NAME flow_ym_za4_num COMMAND compare_ym_za4_num)
  add_test(NAME flow_zq_num COMMAND compare_zq_num)
  add_test(NAME flow_flavour_split COMMAND compare_flavour_split)
  add_test(NAME flow_flavour_ingroup COMMAND compare_flavour_ingroup)
  add_test(NAME flow_gluon_condensate COMMAND compare_gluon_condensate)
  add_test(NAME color_su3 COMMAND sun_contract 3)

  # ---- CODEGEN test (opt-in) -------------------------------------------------------------------
  # Every test above runs a C++ binary against the PRE-GENERATED headers in tests/gen; Wolfram is
  # never invoked, so none of them can see a change to mathematica/. This one closes that gap: it
  # regenerates the kernels from source, rebuilds, and re-runs the flow tests. It needs Wolfram and
  # takes minutes, hence opt-in:
  #
  #     cmake -B build -DNUMTRACER_TEST_CODEGEN=ON && ctest --test-dir build -L codegen
  #
  # It carries the label `codegen` and re-invokes ctest with `-LE codegen`, so it cannot recurse into
  # itself; for the same reason the install-time run below excludes the label.
  option(NUMTRACER_TEST_CODEGEN
         "Add a ctest that regenerates the kernels via Wolfram and re-runs the flow tests (slow)" OFF)
  if(NUMTRACER_TEST_CODEGEN)
    find_program(NT_WOLFRAM_EXE NAMES wolfram wolframscript math)
    if(NOT NT_WOLFRAM_EXE)
      message(WARNING "NUMTRACER_TEST_CODEGEN=ON but no Wolfram kernel found; skipping the codegen test")
    else()
      add_test(NAME codegen_regen
               COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/tests/gen/regen_check.sh" --ctest
               WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
      set_tests_properties(codegen_regen PROPERTIES
        LABELS codegen
        RUN_SERIAL TRUE          # regenerates in the source tree and rebuilds: must not race
        TIMEOUT 5400)
      set_property(TEST codegen_regen PROPERTY ENVIRONMENT "NT_BUILD_DIR=${CMAKE_BINARY_DIR}")
    endif()
  endif()

  # Run the test suite as part of `cmake --install` (and `make/ninja install`, which build
  # `all` first so the test binaries exist). A failing test aborts the install with a fatal
  # error. Off-switch: -DNUMTRACER_TEST_ON_INSTALL=OFF. Tests run against the build-tree
  # binaries, so this verifies the build we are about to install, independent of the prefix.
  option(NUMTRACER_TEST_ON_INSTALL "Run ctest as part of the install target" ON)
  if(NUMTRACER_TEST_ON_INSTALL)
    install(CODE "
      message(STATUS \"NumTracer: running the test suite before finishing install\")
      execute_process(
        COMMAND \"${CMAKE_CTEST_COMMAND}\" --test-dir \"${CMAKE_BINARY_DIR}\"
                --output-on-failure --build-config \"\$<CONFIG>\" -LE codegen
        RESULT_VARIABLE _nt_ctest_rc)
      if(NOT _nt_ctest_rc EQUAL 0)
        message(FATAL_ERROR \"NumTracer: tests failed (ctest exit \${_nt_ctest_rc}); install aborted\")
      endif()")
  endif()
endif()

# ##############################################################################
# GPU integration tests (CUDA + GSL; OFF by default so the core stays dependency-free)
# ##############################################################################
#
# Full 4D loop integrals of the generated ZA3/ZA4 kernels over a grid of external momenta,
# DiFfRG-style two-phase map/reduce on the GPU (tests/gpu/integrator.cuh — an internal test
# tool, not part of the shipped public headers).
# Needs the __host__ __device__ kernel variants (tests/gen/ZA{3,4}_gpu_kernel*.hh, regenerate
# with tests/gen/gen_qcd_za{3,4}_gpu.wls), the CUDA toolkit and GSL (quadrature nodes/weights).
# Configure with -DNUMTRACER_GPU_TESTS=ON; if the default host compiler is too new for nvcc,
# add -DCMAKE_CUDA_HOST_COMPILER=g++-14. clang CUDA also works (-DCMAKE_CUDA_COMPILER=clang++).
option(NUMTRACER_GPU_TESTS "Build the CUDA GPU loop-integral integration tests" OFF)
# Same tests expressed in Kokkos (CUDA backend) instead of raw CUDA + cub, for a framework
# overhead comparison. Uses the DiFfRG-bundled Kokkos install (override NT_KOKKOS_ROOT).
option(NUMTRACER_KOKKOS_TESTS "Build the Kokkos twins of the GPU integration tests" OFF)

if(NUMTRACER_GPU_TESTS OR NUMTRACER_KOKKOS_TESTS)
  find_package(GSL REQUIRED)
  # CPU same-quadrature reference (the committed HOST kernels), shared by the CUDA and Kokkos
  # test mains and always compiled as plain C++ by the configured CXX compiler.
  add_library(nt_gpu_ref STATIC tests/gpu/ref_za3.cpp tests/gpu/ref_za4.cpp)
  # tests/gpu BEFORE tests/gen: the generated kernels' `#include "shim.hpp"` must resolve to
  # the host/device shim, not tests/refshim's host-only one (which is not on this path at all).
  target_include_directories(nt_gpu_ref PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tests/gpu
                                               ${CMAKE_CURRENT_SOURCE_DIR}/tests/gen)
  target_link_libraries(nt_gpu_ref PUBLIC NumTracer::NumTracer GSL::gsl)
endif()

if(NUMTRACER_GPU_TESTS)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)
  endif()
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 20)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)

  function(nt_gpu_test name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE nt_gpu_ref)
    target_compile_options(${name} PRIVATE
      # --expt-relaxed-constexpr lets __device__ code call the constexpr host functions that read
      # the typed-out SU(3) tables (a permission, not a budget). No raised constexpr op/step limits
      # are needed: the SU(N) folding happens in the build-time generator, so the emitted kernel
      # only does table lookups — verified to compile at default limits (nvcc + gcc host pass).
      "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CUDA_COMPILER_ID:NVIDIA>>:--expt-relaxed-constexpr;-O3;-lineinfo>"
      "$<$<AND:$<COMPILE_LANGUAGE:CUDA>,$<CUDA_COMPILER_ID:Clang>>:-O3>")
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES LABELS gpu)
  endfunction()

  nt_gpu_test(integrate_za3_gpu tests/gpu/integrate_za3_gpu.cu)
  nt_gpu_test(integrate_za4_gpu tests/gpu/integrate_za4_gpu.cu)
endif()

if(NUMTRACER_KOKKOS_TESTS)
  set(NT_KOKKOS_ROOT "$ENV{HOME}/.local/share/DiFfRG/bundled" CACHE PATH
      "Kokkos install prefix (default: the DiFfRG bundled install)")
  find_package(Kokkos REQUIRED PATHS ${NT_KOKKOS_ROOT} NO_DEFAULT_PATH)

  function(nt_kokkos_test name src)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE nt_gpu_ref Kokkos::kokkos)
    # Kokkos::kokkos injects --expt-relaxed-constexpr, --extended-lambda and the arch flags; the
    # emitted kernel only reads the typed-out SU(3) tables, so no raised constexpr budget is needed.
    # route this target's C++ compiles through kokkos_launch_compiler -> nvcc_wrapper(g++-14)
    kokkos_compilation(TARGET ${name})
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES LABELS "gpu;kokkos")
  endfunction()

  nt_kokkos_test(integrate_za3_kokkos tests/gpu/integrate_za3_kokkos.cpp)
  nt_kokkos_test(integrate_za4_kokkos tests/gpu/integrate_za4_kokkos.cpp)
endif()

# ##############################################################################
# Documentation
# ##############################################################################
#
# A Sphinx site bundling the C++ API (Doxygen wrapped by doxysphinx) and the
# narrative pages. It is orchestrated by documentation/build.sh, which needs
# Doxygen plus a Python 3.9–3.12 environment with the packages in
# documentation/requirements.txt. OFF by default so an ordinary build never
# needs the docs toolchain; we also degrade gracefully when it is missing.
option(NumTracer_DOCUMENTATION "Build the documentation site" OFF)

if(NumTracer_DOCUMENTATION)
  find_package(Doxygen)
  find_program(SPHINX_BUILD_EXECUTABLE sphinx-build)
  find_program(DOXYSPHINX_EXECUTABLE doxysphinx)

  if(DOXYGEN_FOUND AND SPHINX_BUILD_EXECUTABLE AND DOXYSPHINX_EXECUTABLE)
    message(STATUS "Documentation target configured (Doxygen + doxysphinx + Sphinx)")
    add_custom_target(
      documentation
      COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/documentation/build.sh
              ${CMAKE_CURRENT_BINARY_DIR}/documentation
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/documentation
      COMMENT "Building NumTracer documentation (C++ API + narrative)"
      VERBATIM)
  else()
    # Report exactly which tools are missing so the user knows what to supply.
    set(_doc_missing "")
    if(NOT DOXYGEN_FOUND)
      list(APPEND _doc_missing "doxygen")
    endif()
    if(NOT SPHINX_BUILD_EXECUTABLE)
      list(APPEND _doc_missing "sphinx-build")
    endif()
    if(NOT DOXYSPHINX_EXECUTABLE)
      list(APPEND _doc_missing "doxysphinx")
    endif()
    string(REPLACE ";" ", " _doc_missing_str "${_doc_missing}")
    message(WARNING
      "The 'documentation' target was NOT created: missing ${_doc_missing_str}.\n"
      "  * doxygen — install via your package manager (e.g. 'pacman -S doxygen').\n"
      "  * sphinx-build / doxysphinx — come from the Python toolchain. doxysphinx requires\n"
      "    Python 3.9-3.12 (NOT 3.13+), so create such an environment and install the deps:\n"
      "        python3.12 -m venv ${CMAKE_CURRENT_SOURCE_DIR}/documentation/.venv-docs\n"
      "        source ${CMAKE_CURRENT_SOURCE_DIR}/documentation/.venv-docs/bin/activate\n"
      "        pip install -r ${CMAKE_CURRENT_SOURCE_DIR}/documentation/requirements.txt\n"
      "    then re-run cmake with that environment ACTIVE so the tools are on PATH.\n"
      "  Or skip cmake entirely and build the docs directly:\n"
      "        bash ${CMAKE_CURRENT_SOURCE_DIR}/documentation/build.sh\n"
      "  Configure with -DNumTracer_DOCUMENTATION=OFF to silence this warning.")
  endif()
endif()
