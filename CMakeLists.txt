cmake_minimum_required(VERSION 3.5)
project(libvc)

include(ExternalProject)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
set(CMAKE_CXX_STANDARD 14)

message("Installing dependencies ... ")
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/dep)
execute_process(COMMAND ${CMAKE_COMMAND} ${CMAKE_SOURCE_DIR}/dep WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/dep)
execute_process(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/dep WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/dep RESULT_VARIABLE build_result)
if(build_result)
    message(FATAL_ERROR "Dependency installation failed.")
endif()
message("Done installing dependencies ... ")

add_library(hts STATIC IMPORTED)
set_property(TARGET hts PROPERTY IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/dep/htslib_project-prefix/src/htslib_project/libhts.a)
set_property(TARGET hts APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_BINARY_DIR}/dep/htslib_project-prefix/src/htslib_project/htslib)

list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR}/dep/xzbuf_project-prefix/src/xzbuf_project-build)
find_package(xzbuf)


set(HEADER_FILES
    include/vc.hpp include/reader.hpp include/m3vcf_reader.hpp include/cmf_reader.hpp include/allele_status.hpp include/vcf_reader.hpp include/varint.hpp include/s1r.hpp include/haplotype_vector.hpp include/genotype_vector.hpp include/variant_iterator.hpp include/portable_endian.hpp)

set(SOURCE_FILES
    src/vc.cpp src/reader.cpp src/m3vcf_reader.cpp include/cmf_reader.hpp src/cmf_reader.cpp include/test_class.hpp src/vcf_reader.cpp src/varint.cpp)

add_library(vc ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(vc hts z xzbuf)
target_include_directories(vc PUBLIC include)

add_executable(vc-test src/test.cpp src/test_class.cpp include/test_class.hpp)
target_link_libraries(vc-test vc)

add_executable(bcf2m3vcf src/bcf2m3vcf.cpp)
target_link_libraries(bcf2m3vcf vc)

add_executable(bcf2cvcf src/bcf2cvcf.cpp)
target_link_libraries(bcf2cvcf vc)

add_executable(cvcf-sample-sort src/cmf_sample_sort.cpp)
target_link_libraries(cvcf-sample-sort vc)
