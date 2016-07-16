cmake_minimum_required(VERSION 3.5)
project(libvc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(include)

set(HEADER_FILES
    include/reader.hpp include/m3vcf_reader.hpp)

set(SOURCE_FILES
    src/reader.cpp src/m3vcf_reader.cpp)

add_library(vc ${SOURCE_FILES} ${HEADER_FILES})

add_executable(vc-test src/test.cpp)
target_link_libraries(vc-test vc)