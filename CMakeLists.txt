project(h-tautau)
cmake_minimum_required(VERSION 3.6)

get_filename_component(AnalysisTools_DIR "${PROJECT_SOURCE_DIR}/../AnalysisTools" ABSOLUTE)
include("${AnalysisTools_DIR}/cmake/include/common.cmake")

file(GLOB_RECURSE ANA_TOOLS_SRC "${AnalysisTools_DIR}/*.cpp")
add_library(AnalysisTools STATIC ${ANA_TOOLS_SRC} "${RootDict}")
add_dependencies(AnalysisTools GenerateRootDict)

file(GLOB_RECURSE HTT_UTILITIES_SRC "${CMSSW_BASE_SRC}/HTT-utilities/*.cc")
add_library(HTT-utilities STATIC ${HTT_UTILITIES_SRC})
set_source_files_properties(${HTT_UTILITIES_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE HH_KINFIT2_SRC "${CMSSW_BASE_SRC}/HHKinFit2/*.cpp")
add_library(HHKinFit2 STATIC ${HH_KINFIT2_SRC})
target_include_directories(HHKinFit2 PRIVATE "${CMSSW_BASE_SRC}/HHKinFit2/HHKinFit2/interface")
set_source_files_properties(${HH_KINFIT2_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE SVFIT_SRC "${CMSSW_BASE_SRC}/TauAnalysis/ClassicSVfit/*.cc")
list(FILTER SVFIT_SRC EXCLUDE REGEX ".*test.*")
add_library(SVfit STATIC ${SVFIT_SRC})
set_source_files_properties(${SVFIT_SRC} PROPERTIES COMPILE_FLAGS "-w")

add_library(h-tautau STATIC ${SOURCE_LIST})

foreach(exe_name ${EXE_LIST})
    target_link_libraries("${exe_name}" h-tautau AnalysisTools HTT-utilities SVfit HHKinFit2)
endforeach()

set_target_properties(Print_Graph Print_Graph_2 Print_Graph_3 Print_TreeBranch2D PROPERTIES EXCLUDE_FROM_ALL 1)
