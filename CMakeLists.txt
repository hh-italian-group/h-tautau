project("h-tautau")
cmake_minimum_required(VERSION 2.6)

include("../AnalysisTools/common.cmake")

file(GLOB OTHER_LIST "Run/config/*")
add_custom_target(others SOURCES ${OTHER_LIST})

add_library("h-tautau" OBJECT ${HEADER_LIST})
