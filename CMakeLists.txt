cmake_minimum_required(VERSION 3.5)
project(libvc)

include(ExternalProject)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

ExternalProject_Add(htslib
                    BUILD_IN_SOURCE 1
                    URL https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2
                    URL_MD5 16d78f90b72f29971b042e8da8be6843
                    #GIT_REPOSITORY "https://github.com/samtools/htslib"
                    #GIT_TAG "1.3.1"
                    UPDATE_COMMAND ""
                    CONFIGURE_COMMAND autoconf && ./configure --prefix=<INSTALL_DIR>
                    BUILD_COMMAND make
                    INSTALL_COMMAND make install)

ExternalProject_Get_Property(htslib SOURCE_DIR)
set(HTSLIB_SRC_DIR ${SOURCE_DIR})

include_directories(include ${HTSLIB_SRC_DIR}/htslib)
link_directories(${HTSLIB_SRC_DIR})



set(HEADER_FILES
    include/vc.hpp include/reader.hpp include/m3vcf_reader.hpp include/cmf_reader.hpp include/allele_status.hpp include/vcf_reader.hpp include/varint.hpp)

set(SOURCE_FILES
    src/vc.cpp src/reader.cpp src/m3vcf_reader.cpp include/cmf_reader.hpp src/cmf_reader.cpp include/test_class.hpp src/vcf_reader.cpp src/varint.cpp)

add_library(vc ${SOURCE_FILES} ${HEADER_FILES})
add_dependencies(vc htslib)

add_executable(vc-test src/test.cpp src/test_class.cpp include/test_class.hpp)
target_link_libraries(vc-test vc)
target_link_libraries(vc-test "hts")
target_link_libraries(vc-test "z")

add_executable(bcf2m3vcf src/bcf2m3vcf.cpp)
target_link_libraries(bcf2m3vcf vc)
target_link_libraries(bcf2m3vcf "hts")
target_link_libraries(bcf2m3vcf "z")

add_executable(bcf2cvcf src/bcf2cvcf.cpp)
target_link_libraries(bcf2cvcf vc)
target_link_libraries(bcf2cvcf "hts")
target_link_libraries(bcf2cvcf "z")


add_executable(cvcf-sample-sort src/cmf_sample_sort.cpp)
target_link_libraries(cvcf-sample-sort vc)
target_link_libraries(cvcf-sample-sort "hts")
target_link_libraries(cvcf-sample-sort "z")