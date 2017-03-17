cmake_minimum_required(VERSION 3.2)
project(foo)
message("Dep ROOT: ${CMAKE_ROOT}")
include(ExternalProject)


ExternalProject_Add(htslib_project
                    BUILD_IN_SOURCE 1
                    URL https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2
                    URL_MD5 16d78f90b72f29971b042e8da8be6843
                    #GIT_REPOSITORY "https://github.com/samtools/htslib"
                    #GIT_TAG "1.3.1"
                    UPDATE_COMMAND ""
                    CONFIGURE_COMMAND "" #autoconf && ./configure --prefix=<INSTALL_DIR>
                    BUILD_COMMAND make
                    INSTALL_COMMAND "")


ExternalProject_Add(xzbuf_project
                    GIT_REPOSITORY "https://github.com/jonathonl/xzbuf"
                    GIT_TAG "develop"
                    INSTALL_COMMAND "")