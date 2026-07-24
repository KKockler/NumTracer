# NumTracerNumtrace.cmake — offline kernel generation as a build target.
#
#   numtracer_add_numtrace(<flows_target> <flow_dir>)
#
# Emitting a flow from Wolfram writes three things: the generator sources (into a gen/ directory that
# is normally gitignored), a probe source when the flow is complex, and a per-flow
# <flow_dir>/<Name>/numtrace.json manifest. This function turns the manifests into a `numtrace` target
# that compiles and runs whatever still needs generating, and makes <flows_target> depend on it — so
# the flows library can never be built against stale or missing kernels.
#
# The switch is the manifest's "generated" field, not a file timestamp. It has to be: kernels.hh is
# committed while gen/ is not, so a fresh clone has the outputs and none of the inputs, and a
# timestamp rule would either demand sources that do not exist or silently accept stale ones. A flow
# at 1 contributes nothing here — if every flow is at 1, `numtrace` is an empty target and does
# nothing, which is the normal state of a checked-out tree.
#
# Parallelism comes from CMake itself: each flow's generator is an ordinary executable target, so
# `make -j8 numtrace` spreads eight compiles across all flows at once rather than running one
# per-flow xargs after another, which is what the Wolfram path did.
#
# Per flow the chain is
#
#   gen_<ns>_num.cpp + _u*.cpp  --compile+link-->  numtrace_gen_<ns>
#                                      --run-->   <Name>/kernels.hh
#                                   --compile-->  numtrace_probe_<ns>      (complex flows only)
#                                      --run-->   <Name>/numtrace_verdict.hh
#
# and the last step in the chain also flips the switch, so a flow is marked done exactly once and
# only after everything it owes has succeeded. The verdict header carries the imaginary-part probe's
# three-way answer; kernel.hh holds all three bodies under an #if and the preprocessor picks one.

if(COMMAND numtracer_add_numtrace)
  return()
endif()

function(numtracer_add_numtrace flows_target flow_dir)
  set(_driver "${NumTracer_NUMTRACE_DIR}/NumTracerNumtraceRun.cmake")
  if(NOT EXISTS "${_driver}")
    message(FATAL_ERROR "numtrace: driver script not found at ${_driver}")
  endif()

  file(GLOB _manifests "${flow_dir}/*/numtrace.json")

  # A re-emit from Wolfram rewrites the manifests (switch back to 0). Registering them here makes the
  # next build reconfigure by itself, so the new work is picked up without anyone remembering to
  # re-run cmake. The glob is registered too, so an entirely new flow directory also triggers it.
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_manifests})

  set(_outputs "")
  set(_pending "")

  foreach(_manifest IN LISTS _manifests)
    get_filename_component(_flowdir "${_manifest}" DIRECTORY)
    file(READ "${_manifest}" _json)

    string(JSON _generated GET "${_json}" "generated")
    if(_generated)
      continue()
    endif()

    string(JSON _name      GET "${_json}" "name")
    string(JSON _ns        GET "${_json}" "namespace")
    string(JSON _gendir    GET "${_json}" "gen_dir")
    string(JSON _generator GET "${_json}" "generator")
    string(JSON _decorator GET "${_json}" "decorator")
    string(JSON _mainopt   GET "${_json}" "main_opt")
    string(JSON _fullpar   GET "${_json}" "full_parallel")
    string(JSON _complex   GET "${_json}" "complex")
    string(JSON _kernels   GET "${_json}" "kernels")

    # Manifest paths are normally relative to the flow directory, but tolerate absolute ones: a path
    # that could not be made relative at emit time (e.g. a not-yet-created file whose symlinked parent
    # resolved differently) is written out absolute, and must not be re-prefixed.
    if(NOT IS_ABSOLUTE "${_gendir}")
      set(_gendir "${_flowdir}/${_gendir}")
    endif()
    get_filename_component(_gendir "${_gendir}" ABSOLUTE)

    # unit sources
    string(JSON _nunits LENGTH "${_json}" "units")
    set(_sources "${_gendir}/${_generator}")
    math(EXPR _last "${_nunits} - 1")
    if(_nunits GREATER 0)
      foreach(_i RANGE ${_last})
        string(JSON _unit GET "${_json}" "units" ${_i})
        list(APPEND _sources "${_gendir}/${_unit}")
      endforeach()
    endif()

    # A committed switch of 0 with no sources on disk is unbuildable, and the resulting CMake error
    # would otherwise be an obscure "cannot find source file". Say what actually has to happen.
    foreach(_src IN LISTS _sources)
      if(NOT EXISTS "${_src}")
        message(FATAL_ERROR
          "numtrace: flow '${_name}' is marked ungenerated but its generator source is missing:\n"
          "  ${_src}\n"
          "Re-run the flow notebook to re-emit it, or set \"generated\": 1 in\n"
          "  ${_manifest}\n"
          "if the committed kernels are known to be current.")
      endif()
    endforeach()

    # ---- the generator -------------------------------------------------------------------------
    set(_gentgt "numtrace_gen_${_ns}")
    add_executable(${_gentgt} EXCLUDE_FROM_ALL ${_sources})
    target_link_libraries(${_gentgt} PRIVATE NumTracer::NumTracer)
    target_compile_options(${_gentgt} PRIVATE -ftemplate-depth=4000)
    set_target_properties(${_gentgt} PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

    # The main TU carries the contraction loop and must stay optimised or the RUN slows ~10x; the
    # net-builder units are pure straight-line codegen where -O0 is both faster to compile and
    # runtime-neutral. -fno-exceptions -fno-rtti on the units drops the exception-cleanup landing pads
    # for their thousands of destructible braced-init temporaries (~90% of the -O0 unit compile); the
    # unit code never throws. Same split and flags the Wolfram generator compile uses.
    set_source_files_properties("${_gendir}/${_generator}"
      TARGET_DIRECTORY ${_gentgt} PROPERTIES COMPILE_OPTIONS "${_mainopt}")
    foreach(_src IN LISTS _sources)
      if(NOT _src STREQUAL "${_gendir}/${_generator}")
        set_source_files_properties("${_src}"
          TARGET_DIRECTORY ${_gentgt} PROPERTIES COMPILE_OPTIONS "-O0;-fno-exceptions;-fno-rtti")
      endif()
    endforeach()

    if(IS_ABSOLUTE "${_kernels}")
      set(_kernels_hh "${_kernels}")
    else()
      set(_kernels_hh "${_flowdir}/${_kernels}")
    endif()
    set(_run_cmd
      COMMAND ${CMAKE_COMMAND}
        -DMODE=run
        "-DGEN=$<TARGET_FILE:${_gentgt}>"
        "-DOUT=${_kernels_hh}"
        "-DNS=${_ns}"
        "-DDECOR=${_decorator}"
        "-DFULLPAR=${_fullpar}"
        -P "${_driver}")

    if(_complex)
      # kernels.hh first, then the probe that reads it; the switch flips at the end of the chain.
      add_custom_command(
        OUTPUT "${_kernels_hh}"
        ${_run_cmd}
        DEPENDS ${_gentgt}
        COMMENT "NumTracer: tracing ${_name}"
        VERBATIM)

      string(JSON _probe   GET "${_json}" "probe")
      string(JSON _macro   GET "${_json}" "verdict_macro")
      string(JSON _verdict GET "${_json}" "verdict")

      set(_probetgt "numtrace_probe_${_ns}")
      add_executable(${_probetgt} EXCLUDE_FROM_ALL "${_gendir}/${_probe}")
      target_include_directories(${_probetgt} PRIVATE "${_flowdir}")
      target_compile_options(${_probetgt} PRIVATE -O1 -w)
      set_target_properties(${_probetgt} PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
      # the probe #includes the traces header, so it cannot be compiled before the generator has run
      set_source_files_properties("${_gendir}/${_probe}"
        TARGET_DIRECTORY ${_probetgt} PROPERTIES OBJECT_DEPENDS "${_kernels_hh}")

      if(IS_ABSOLUTE "${_verdict}")
        set(_verdict_hh "${_verdict}")
      else()
        set(_verdict_hh "${_flowdir}/${_verdict}")
      endif()
      add_custom_command(
        OUTPUT "${_verdict_hh}"
        COMMAND ${CMAKE_COMMAND}
          -DMODE=probe
          "-DPROBE=$<TARGET_FILE:${_probetgt}>"
          "-DOUT=${_verdict_hh}"
          "-DMACRO=${_macro}"
          "-DFLOW=${_name}"
          -P "${_driver}"
        DEPENDS ${_probetgt} "${_kernels_hh}"
        COMMENT "NumTracer: probing ${_name}"
        VERBATIM)

      set(_last "${_verdict_hh}")
      list(APPEND _outputs "${_kernels_hh}" "${_verdict_hh}")
    else()
      add_custom_command(
        OUTPUT "${_kernels_hh}"
        ${_run_cmd}
        DEPENDS ${_gentgt}
        COMMENT "NumTracer: tracing ${_name}"
        VERBATIM)
      set(_last "${_kernels_hh}")
      list(APPEND _outputs "${_kernels_hh}")
    endif()

    # The switch flip gets its OWN output, a stamp in the build tree, rather than riding along as a
    # trailing COMMAND of the rule that produces kernels.hh. Make deletes a rule's output when the
    # recipe fails, so a failure in the flip would otherwise delete the committed traces header that
    # the run step had already written. The stamp is expendable; the header is not.
    set(_stampdir "${CMAKE_CURRENT_BINARY_DIR}/numtrace")
    set(_stamp "${_stampdir}/${_ns}.marked")
    file(MAKE_DIRECTORY "${_stampdir}")
    add_custom_command(
      OUTPUT "${_stamp}"
      COMMAND ${CMAKE_COMMAND} -DMODE=mark "-DMANIFEST=${_manifest}" -P "${_driver}"
      COMMAND ${CMAKE_COMMAND} -E touch "${_stamp}"
      DEPENDS "${_last}"
      COMMENT "NumTracer: ${_name} done"
      VERBATIM)
    list(APPEND _outputs "${_stamp}")

    list(APPEND _pending "${_name}")
  endforeach()

  # `numtrace` is always defined, so `make numtrace` and the dependency below work whether or not
  # anything is pending. With every switch at 1 it is an empty target. A project with more than one
  # flow directory calls this function once per directory, so the second call hangs its work off a
  # per-directory sub-target rather than redefining `numtrace`.
  if(TARGET numtrace)
    string(MAKE_C_IDENTIFIER "numtrace_${flows_target}" _sub)
    add_custom_target(${_sub} DEPENDS ${_outputs})
    add_dependencies(numtrace ${_sub})
  else()
    add_custom_target(numtrace DEPENDS ${_outputs})
  endif()
  add_dependencies(${flows_target} numtrace)

  if(_pending)
    list(LENGTH _pending _n)
    string(REPLACE ";" " " _pretty "${_pending}")
    message(STATUS "numtrace: ${_n} flow(s) to generate: ${_pretty}")
  elseif(_manifests)
    message(STATUS "numtrace: all flows already generated")
  else()
    # No manifests at all: either the flows predate offline generation (they were emitted with
    # "Offline" -> False, so their kernels are already committed) or nothing has been emitted yet.
    # Either way there is nothing for this target to do.
    message(STATUS "numtrace: no numtrace manifests under ${flow_dir} — nothing to generate")
  endif()
endfunction()
