cmake_minimum_required(VERSION 3.2)
project(htslib VERSION 1.3.1)

#execute_process(COMMAND ./configure --disable-libcurl --disable-lzma --disable-bz2 --prefix=${CMAKE_INSTALL_PREFIX} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

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

#install(DIRECTORY htslib DESTINATION include)

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
