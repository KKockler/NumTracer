set(QCDFlows_SOURCES 
    ${CMAKE_CURRENT_SOURCE_DIR}/ZA/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZA/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZA/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZA3/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZA3/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZA3/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZA4/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZA4/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZA4/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAAqbq1/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAAqbq1/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAAqbq1/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAAqbq2/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAAqbq2/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAAqbq2/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAcbc/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAcbc/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAcbc/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAqbq1/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAqbq1/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAqbq1/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAqbq4/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAqbq4/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAqbq4/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAqbq7/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAqbq7/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ZAqbq7/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Zc/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Zc/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Zc/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Zq/src/constructor.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Zq/src/CT_get.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/Zq/src/CT_map_1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/flows.cc)

find_package(NumTracer REQUIRED HINTS ~/.local/share/NumTracer)

add_library(QCDFlows STATIC ${QCDFlows_SOURCES})
add_library(QCDFlows_nowarn INTERFACE)
target_compile_options(QCDFlows_nowarn INTERFACE -Wno-unused-parameter -Wno-unused-variable)
target_link_libraries(QCDFlows DiFfRG::DiFfRG QCDFlows_nowarn NumTracer::NumTracer)

numtracer_add_numtrace(QCDFlows ${CMAKE_CURRENT_SOURCE_DIR})

set_target_properties(
  QCDFlows
  PROPERTIES
    UNITY_BUILD OFF
    UNITY_BUILD_MODE GROUP)

foreach(src IN LISTS QCDFlows_SOURCES)
  if(src MATCHES "/src/CT_map_[0-9]+\.cc$")
    set_source_files_properties("${src}" PROPERTIES UNITY_GROUP "all_ct_map")
  else()
    set_source_files_properties("${src}" PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
  endif()
endforeach()