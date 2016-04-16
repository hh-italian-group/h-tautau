project("h-tautau")
cmake_minimum_required(VERSION 2.6)

set(AnalysisTools_DIR "${PROJECT_SOURCE_DIR}/../AnalysisTools")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${AnalysisTools_DIR}/cmake/modules")
include("${AnalysisTools_DIR}/common.cmake")

file(GLOB OTHER_LIST "Run/config/*")
add_custom_target(others SOURCES ${OTHER_LIST})

add_library("h-tautau" OBJECT ${HEADER_LIST})
