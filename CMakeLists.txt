project(h-tautau)
cmake_minimum_required(VERSION 3.0)

get_filename_component(AnalysisTools_DIR "${PROJECT_SOURCE_DIR}/../AnalysisTools" ABSOLUTE)
include("${AnalysisTools_DIR}/cmake/include/common.cmake")

file(GLOB_RECURSE HTT_UTILITIES_SRC "${CMSSW_BASE_SRC}/HTT-utilities/*.cc")
add_library(HTT-utilities STATIC ${HTT_UTILITIES_SRC})
set_source_files_properties(${HTT_UTILITIES_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE SVFIT_SRC "${CMSSW_BASE_SRC}/TauAnalysis/SVfitStandalone/*.cc")
add_library(SVfit STATIC ${SVFIT_SRC})
set_source_files_properties(${SVFIT_SRC} PROPERTIES COMPILE_FLAGS "-w")

foreach(exe_name ${EXE_LIST})
    target_link_libraries("${exe_name}" HTT-utilities SVfit)
endforeach()

set_target_properties(Print_Graph Print_Graph_2 Print_Graph_3 Print_Stack Print_TreeBranch2D PROPERTIES EXCLUDE_FROM_ALL 1)

add_library("h-tautau" OBJECT ${HEADER_LIST} ${SOURCE_LIST} ${EXE_SOURCE_LIST} ${SCRIPT_LIST} ${CONFIG_LIST})
set_target_properties("h-tautau" PROPERTIES EXCLUDE_FROM_ALL 1)
