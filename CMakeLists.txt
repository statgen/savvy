cmake_minimum_required(VERSION 3.5)
project(libsavvy VERSION 1.0.0)
include(CMakePackageConfigHelpers)
include(ExternalProject)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
set(CMAKE_CXX_STANDARD 14)

message("Installing dependencies ... ")
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/dep)
execute_process(COMMAND ${CMAKE_COMMAND} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX} ${CMAKE_SOURCE_DIR}/dep WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/dep)
execute_process(COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/dep WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/dep RESULT_VARIABLE build_result)
if(build_result)
    message(FATAL_ERROR "Dependency installation failed.")
endif()
message("Done installing dependencies ... ")

add_library(hts STATIC IMPORTED)
set_property(TARGET hts PROPERTY IMPORTED_LOCATION ${CMAKE_INSTALL_PREFIX}/lib/libhts.a)
set_property(TARGET hts APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_INSTALL_PREFIX}/include/htslib)

#list(APPEND CMAKE_PREFIX_PATH ${CMAKE_INSTALL_PREFIX})
find_package(xzbuf REQUIRED)


set(HEADER_FILES
    include/savvy/savvy.hpp include/savvy/reader.hpp include/savvy/m3vcf_reader.hpp include/savvy/sav_reader.hpp include/savvy/allele_status.hpp include/savvy/vcf_reader.hpp include/savvy/varint.hpp include/savvy/s1r.hpp include/savvy/allele_vector.hpp include/savvy/genotype_vector.hpp include/savvy/variant_iterator.hpp include/savvy/portable_endian.hpp include/savvy/region.hpp include/savvy/compressed_vector.hpp)

set(SOURCE_FILES
    src/savvy.cpp src/reader.cpp src/m3vcf_reader.cpp include/savvy/sav_reader.hpp src/sav_reader.cpp include/savvy/test_class.hpp src/vcf_reader.cpp src/varint.cpp)

add_library(savvy ${SOURCE_FILES} ${HEADER_FILES})
target_link_libraries(savvy hts z xzbuf)
target_include_directories(savvy PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)

add_executable(savvy-test src/test.cpp src/test_class.cpp include/savvy/test_class.hpp)
target_link_libraries(savvy-test savvy)

add_executable(bcf2m3vcf src/bcf2m3vcf.cpp)
target_link_libraries(bcf2m3vcf savvy)

add_executable(bcf2sav src/bcf2sav.cpp)
target_link_libraries(bcf2sav savvy)

add_executable(sav-sample-sort src/sav_sample_sort.cpp)
target_link_libraries(sav-sample-sort savvy)

install(DIRECTORY include/savvy DESTINATION include)
install(TARGETS savvy EXPORT ${PROJECT_NAME}-config
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

install(EXPORT ${PROJECT_NAME}-config DESTINATION share/${PROJECT_NAME})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake DESTINATION share/${PROJECT_NAME})

export(EXPORT ${PROJECT_NAME}-config)