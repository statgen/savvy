cmake_minimum_required(VERSION 3.2)
project(htslib VERSION 1.3.1)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I${dir}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${dir}")
endforeach()
message("dirs: ${dirs}")
message("cxxflags: ${CMAKE_CXX_FLAGS}")

#list(APPEND CMAKE_C_FLAGS ${dirs})
#list(APPEND CMAKE_CXX_FLAGS ${dirs})

#execute_process(COMMAND ./configure --disable-libcurl --disable-lzma --disable-bz2 --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set(ENV{CFLAGS}  "${CMAKE_C_FLAGS}")
set(ENV{CXXFLAGS} "${CMAKE_CXX_FLAGS}")
set(ENV{LDFLAGS} "-L${CMAKE_PREFIX_PATH}/lib")

#set(ENV{CFLAGS} "${CMAKE_C_FLAGS} $ENV{CFLAGS}")
#set(ENV{CXXFLAGS} "${CMAKE_CXX_FLAGS} $ENV{CXXFLAGS}")
#set(ENV{AM_CFLAGS} -I${CGET_PREFIX}/include)
#set(ENV{CPPFLAGS} -I${CGET_PREFIX}/include)

execute_process(COMMAND ./configure --disable-libcurl --disable-lzma --disable-bz2 --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_custom_command(OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hts${CMAKE_SHARED_LIBRARY_SUFFIX}" "${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}hts${CMAKE_STATIC_LIBRARY_SUFFIX}"
                   COMMAND $(MAKE)
                   WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                   COMMENT "Building htslib ...")
add_custom_target(hts ALL DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}hts${CMAKE_SHARED_LIBRARY_SUFFIX}" "${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}hts${CMAKE_STATIC_LIBRARY_SUFFIX}")

#add_custom_target(hts ALL
                  #COMMAND $(MAKE)
                  #WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" 
                  #COMMENT "Building htslib ...")

install(DIRECTORY htslib DESTINATION include)

if (BUILD_SHARED_LIBS)
    install(FILES
            ${CMAKE_SHARED_LIBRARY_PREFIX}hts${CMAKE_SHARED_LIBRARY_SUFFIX}
            DESTINATION lib)
    install(FILES
            ${CMAKE_SHARED_LIBRARY_PREFIX}hts.1${CMAKE_SHARED_LIBRARY_SUFFIX}
            ${CMAKE_SHARED_LIBRARY_PREFIX}hts${CMAKE_SHARED_LIBRARY_SUFFIX}.1
            DESTINATION lib
            OPTIONAL)
else()
    install(FILES
            ${CMAKE_STATIC_LIBRARY_PREFIX}hts${CMAKE_STATIC_LIBRARY_SUFFIX}
            DESTINATION lib)
endif()
