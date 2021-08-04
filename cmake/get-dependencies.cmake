#!/usr/bin/env cmake -P

set(ARGS)
foreach(i RANGE 4 ${CMAKE_ARGC})
    list(APPEND ARGS ${CMAKE_ARGV${i}})
endforeach()

set(_PREFIX ${CMAKE_ARGV3})

file(DOWNLOAD https://raw.githubusercontent.com/pfultz2/cmake-get/03ef767603d21464fd7e89ea8cb6007a90c472d2/share/cmake/cmakeget/CMakeGet.cmake ${CMAKE_CURRENT_LIST_DIR}/CMakeGet.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/CMakeGet.cmake)

get_filename_component(PREFIX ${_PREFIX} ABSOLUTE)

cmake_get_from(requirements.txt PREFIX ${PREFIX} CMAKE_ARGS ${ARGS})
