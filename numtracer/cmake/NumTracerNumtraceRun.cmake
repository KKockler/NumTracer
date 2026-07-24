# NumTracerNumtraceRun.cmake — build-step driver for the `numtrace` target (cmake -P).
#
# Two modes, selected with -DMODE=:
#
#   run   run a flow's generator binary and commit its stdout as the straight-line traces header.
#         -DGEN=<binary> -DOUT=<kernels.hh> -DNS=<namespace> -DDECOR=<decorator> [-DFULLPAR=ON]
#
#   probe run a flow's imaginary-part probe, which writes the verdict header, and report the verdict
#         in words. The probe's raw measurement line is captured rather than dumped into the build log.
#         -DPROBE=<binary> -DOUT=<numtrace_verdict.hh> -DMACRO=<name> -DFLOW=<name>
#
#   mark  flip a flow's numtrace.json switch to "generated": 1.
#         -DMANIFEST=<numtrace.json>
#
# `run` mirrors the safety checks the Wolfram path used to do inline: the generator writes to a temp
# file, which is validated (exit code 0 AND non-trivial size) before replacing the committed header —
# a crashed generator must never truncate a header that is checked into git.

cmake_minimum_required(VERSION 3.19)

if(MODE STREQUAL "run")
  foreach(v GEN OUT NS DECOR)
    if(NOT DEFINED ${v})
      message(FATAL_ERROR "numtrace run: -D${v}= is required")
    endif()
  endforeach()

  set(_tmp "${OUT}.tmp")
  set(_args -n "${NS}" -d "${DECOR}")
  if(FULLPAR)
    # heavy nets are reduced/rebased concurrently: faster codegen, higher peak RAM
    list(APPEND _args -p)
  endif()

  execute_process(
    COMMAND "${GEN}" ${_args}
    OUTPUT_FILE "${_tmp}"
    ERROR_VARIABLE _err
    RESULT_VARIABLE _rc)

  set(_sz 0)
  if(EXISTS "${_tmp}")
    file(SIZE "${_tmp}" _sz)
  endif()

  if(NOT _rc EQUAL 0 OR _sz LESS 64)
    file(REMOVE "${_tmp}")
    message(FATAL_ERROR
      "numtrace: generator for '${NS}' failed (exit ${_rc}, ${_sz} bytes).\n"
      "The committed ${OUT} was left intact.\n${_err}")
  endif()

  file(RENAME "${_tmp}" "${OUT}")
  message(STATUS "NumTracer: wrote ${OUT} (${_sz} bytes)")

elseif(MODE STREQUAL "probe")
  foreach(v PROBE OUT MACRO FLOW)
    if(NOT DEFINED ${v})
      message(FATAL_ERROR "numtrace probe: -D${v}= is required")
    endif()
  endforeach()

  execute_process(
    COMMAND "${PROBE}" -o "${OUT}" -m "${MACRO}"
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _err
    RESULT_VARIABLE _rc)

  # Any failure is fatal: the verdict picks a preprocessor branch in a committed header, so a
  # probe that did not actually run must never be mistaken for "it came out complex".
  if(NOT _rc EQUAL 0 OR NOT EXISTS "${OUT}")
    message(FATAL_ERROR
      "numtrace: imaginary-part probe for '${FLOW}' failed (exit ${_rc}).\n${_err}${_out}")
  endif()

  # stdout is "mim mdiff mre mrim mrdiff npoints verdict"; report the last two in words and keep the
  # relative imaginary part, which is the number worth seeing when a flow unexpectedly stays complex.
  string(STRIP "${_out}" _out)
  string(REPLACE " " ";" _fields "${_out}")
  list(LENGTH _fields _nf)
  if(_nf EQUAL 7)
    list(GET _fields 3 _relim)
    list(GET _fields 5 _npts)
    list(GET _fields 6 _v)
    if(_v EQUAL 2)
      set(_word "Pure (imaginary coefficients dropped)")
    elseif(_v EQUAL 1)
      set(_word "RePart (re/im split)")
    else()
      set(_word "complex")
    endif()
    message(STATUS "NumTracer: ${FLOW} verdict = ${_word}, rel|Im|=${_relim} over ${_npts} pts")
  else()
    message(STATUS "NumTracer: ${FLOW} probed")
  endif()

elseif(MODE STREQUAL "mark")
  if(NOT DEFINED MANIFEST)
    message(FATAL_ERROR "numtrace mark: -DMANIFEST= is required")
  endif()
  file(READ "${MANIFEST}" _json)
  string(JSON _json SET "${_json}" "generated" 1)
  file(WRITE "${MANIFEST}" "${_json}")
  message(STATUS "NumTracer: ${MANIFEST} -> generated: 1")

else()
  message(FATAL_ERROR "numtrace: unknown -DMODE=${MODE} (expected 'run', 'probe' or 'mark')")
endif()
