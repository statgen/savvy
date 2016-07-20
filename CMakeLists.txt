cmake_minimum_required(VERSION 3.5)
project(libvc)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(include)

set(HEADER_FILES
    include/reader.hpp include/m3vcf_reader.hpp include/cvcf_reader.hpp include/allele_status.hpp)

set(SOURCE_FILES
    src/reader.cpp src/m3vcf_reader.cpp include/cvcf_reader.hpp src/cvcf_reader.cpp include/test_class.hpp)

add_library(vc ${SOURCE_FILES} ${HEADER_FILES})

add_executable(vc-test src/test.cpp src/test_class.cpp include/test_class.hpp)
target_link_libraries(vc-test vc)