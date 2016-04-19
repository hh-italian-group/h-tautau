project(h-tautau)
cmake_minimum_required(VERSION 2.8)

get_filename_component(AnalysisTools_DIR "${PROJECT_SOURCE_DIR}/../AnalysisTools" ABSOLUTE)
include("${AnalysisTools_DIR}/cmake/include/common.cmake")

add_library("h-tautau" OBJECT ${HEADER_LIST})

set_target_properties(Print_Graph Print_Graph_2 Print_Graph_3 Print_Stack Print_SyncPlots Print_TreeBranch2D PROPERTIES EXCLUDE_FROM_ALL 1)
