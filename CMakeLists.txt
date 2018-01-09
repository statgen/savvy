cmake_minimum_required(VERSION 3.2)
project(savvy VERSION 1.0.0)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 11)

#add_library(hts STATIC IMPORTED)
#set_property(TARGET hts PROPERTY IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib/libhts.a)
#set_property(TARGET hts APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_INSTALL_PREFIX}/include)

if(CGET_PREFIX)
    link_directories(${CGET_PREFIX}/lib)
endif()

#list(APPEND CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX})
#find_package(shrinkwrap CONFIG REQUIRED)
find_library(HTS_LIBRARY hts)
find_library(ZLIB_LIBRARY z)
find_library(ZSTD_LIBRARY zstd)
find_package(Threads)

#get_target_property(SHRINKWRAP_LIBS shrinkwrap INTERFACE_LINK_LIBRARIES)

add_definitions(-DSAVVY_VERSION="${PROJECT_VERSION}-beta")

set(HEADER_FILES
    include/savvy/savvy.hpp include/savvy/reader.hpp include/savvy/m3vcf_reader.hpp include/savvy/sav_reader.hpp include/savvy/allele_status.hpp include/savvy/vcf_reader.hpp include/savvy/varint.hpp include/savvy/s1r.hpp include/savvy/site_info.hpp include/savvy/variant_iterator.hpp include/savvy/portable_endian.hpp include/savvy/region.hpp include/savvy/compressed_vector.hpp include/savvy/eigen3_vector.hpp include/savvy/ublas_vector.hpp include/savvy/armadillo_vector.hpp include/savvy/utility.hpp include/savvy/data_format.hpp)

set(SOURCE_FILES
    src/savvy/savvy.cpp src/savvy/reader.cpp src/savvy/m3vcf_reader.cpp include/savvy/sav_reader.hpp src/savvy/sav_reader.cpp include/test/test_class.hpp src/savvy/vcf_reader.cpp src/savvy/varint.cpp src/savvy/utility.cpp src/savvy/site_info.cpp src/savvy/region.cpp)

add_library(savvy ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(savvy ${HTS_LIBRARY} ${ZLIB_LIBRARY} ${ZSTD_LIBRARY} ${CMAKE_THREAD_LIBS_INIT})
target_include_directories(savvy PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

add_executable(sav
               src/sav/main.cpp
               src/sav/export.cpp include/sav/export.hpp
               src/sav/import.cpp include/sav/import.hpp
               src/sav/index.cpp include/sav/index.hpp
               src/sav/merge.cpp include/sav/merge.hpp
               src/sav/sort.cpp include/sav/sort.hpp
               src/sav/utility.cpp include/sav/utility.hpp)
target_link_libraries(sav savvy)

#add_executable(bcf2m3vcf src/sav/bcf2m3vcf.cpp)
#target_link_libraries(bcf2m3vcf savvy)

#add_executable(sav-sample-sort src/sav/sav_sample_sort.cpp)
#target_link_libraries(sav-sample-sort savvy)

#add_executable(savvy-speed-test src/test/savvy_speed_test.cpp)
#target_link_libraries(savvy-speed-test savvy)


if(BUILD_TESTS)
    enable_testing()

    add_definitions(-DSAVVYT_VCF_FILE=\"${CMAKE_CURRENT_SOURCE_DIR}/test_file.vcf\"
                    -DSAVVYT_SAV_FILE_HARD=\"${CMAKE_CURRENT_SOURCE_DIR}/test_file_hard.sav\"
                    -DSAVVYT_SAV_FILE_DOSE=\"${CMAKE_CURRENT_SOURCE_DIR}/test_file_dose.sav\"
                    -DSAVVYT_MARKER_COUNT_HARD=28
                    -DSAVVYT_MARKER_COUNT_DOSE=20)

    add_executable(savvy-test src/test/main.cpp src/test/test_class.cpp include/test/test_class.hpp)
    target_link_libraries(savvy-test savvy)

    add_test(convert_file_test savvy-test convert-file)
    add_test(subset_test savvy-test subset)
    add_test(varint_test savvy-test varint)
endif()

if(BUILD_SLR_EXAMPLES)
    add_executable(slr-examples src/test/slr_examples.cpp)
    target_link_libraries(slr-examples savvy armadillo)

    find_library(OPENCL_LIB OpenCL)

    add_executable(linreg-ttest src/test/linreg_ttest.cpp)
    target_link_libraries(linreg-ttest savvy ${OPENCL_LIB})
endif()

install(DIRECTORY include/savvy DESTINATION include)
install(TARGETS savvy EXPORT ${PROJECT_NAME}-config
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(TARGETS sav RUNTIME DESTINATION bin)

install(EXPORT ${PROJECT_NAME}-config DESTINATION share/${PROJECT_NAME})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake DESTINATION share/${PROJECT_NAME})

export(EXPORT ${PROJECT_NAME}-config)