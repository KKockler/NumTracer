# gpubench — GPU measurement harness for NumTracer-emitted production kernels.
#
# Builds ONE flow tree (a vendored copy, or a rewritten variant of one) plus the timing driver.
# Nothing here reaches back into the QCD codes: NTB_FLOWS points at a tree under kernels/ or
# variants/, so an experiment can never disturb the sources it was copied from.
#
#   cmake -S . -B build/baseline -DNTB_FLOWS=kernels/ym
#   cmake --build build/baseline -j
#   NTB_CONFIG=kernels/ym/parameter.toml ./build/baseline/ntb_bench
#
# Both flow sets build. Nf2's generated headers take DiFfRG::JSONValue, which is a deprecated alias
# for ConfigTree (config_tree.hh:243), so they compile against the current install unchanged.
cmake_minimum_required(VERSION 3.24)
project(gpubench LANGUAGES CXX)

set(NTB_FLOWS "kernels/ym" CACHE STRING "flow tree to build (relative to this dir, or absolute)")

if(NOT IS_ABSOLUTE "${NTB_FLOWS}")
  set(NTB_FLOWS "${CMAKE_CURRENT_SOURCE_DIR}/${NTB_FLOWS}")
endif()
if(NOT EXISTS "${NTB_FLOWS}/flows/flows.hh")
  message(FATAL_ERROR "NTB_FLOWS=${NTB_FLOWS} has no flows/flows.hh")
endif()
message(STATUS "gpubench: flow tree = ${NTB_FLOWS}")

find_package(DiFfRG REQUIRED HINTS ~/.local/share/DiFfRG)
find_package(NumTracer REQUIRED HINTS ~/.local/share/NumTracer)

file(GLOB NTB_FLOW_SOURCES "${NTB_FLOWS}/flows/*/src/*.cc" "${NTB_FLOWS}/flows/flows.cc")

add_library(ntb_flows STATIC ${NTB_FLOW_SOURCES})
# The generated kernels are full of intentionally unused parameters (the [[maybe_unused]] fill
# signature is uniform across flows); the production build silences the same set.
target_compile_options(ntb_flows PRIVATE -Wno-unused-parameter -Wno-unused-variable
                                         -Wno-unused-local-typedefs -Wno-deprecated-declarations)
target_link_libraries(ntb_flows PUBLIC DiFfRG::DiFfRG NumTracer::NumTracer)
target_include_directories(ntb_flows PUBLIC "${NTB_FLOWS}")

# Deliberately NOT calling numtracer_add_numtrace(): the harness measures the committed kernels and
# must never trigger a regeneration as a build side effect.

# Pick the driver by which flow set was vendored. Zq is unique to QCD_Nf2/no_mesons; the two codes
# have different flow classes, interpolator sets and config formats, so they need separate drivers.
if(EXISTS "${NTB_FLOWS}/flows/Zq")
  set(NTB_DRIVER harness/bench_nf2.cc)
else()
  set(NTB_DRIVER harness/bench.cc)
endif()
message(STATUS "gpubench: driver = ${NTB_DRIVER}")

add_executable(ntb_bench ${NTB_DRIVER})
target_link_libraries(ntb_bench PRIVATE ntb_flows)
target_include_directories(ntb_bench PRIVATE "${NTB_FLOWS}")
# JSONValue is a deprecated alias for ConfigTree; the Nf2 headers still use the old spelling.
target_compile_options(ntb_bench PRIVATE -Wno-unused-parameter -Wno-unused-variable
                                         -Wno-deprecated-declarations)
