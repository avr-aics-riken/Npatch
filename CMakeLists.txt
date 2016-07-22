###
###################################################################################
#
# Npatch - Nagata Patch Library
#
# Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
###################################################################################
###


cmake_minimum_required(VERSION 3.0)

enable_language(CXX)

project(Npatch)
# Npatch  install dir.
message(CMAKE_INSTALL_PREFIX: "${CMAKE_INSTALL_PREFIX}")
set(NPT_DIR "${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}")
message("Npatch install dir : ${NPT_DIR}")


# Version info.
set(VERSION "1.0.1")
set( ${PROJECT_NAME}_MAJOR_VERSION 1)
set( ${PROJECT_NAME}_MINOR_VERSION 0)
set( ${PROJECT_NAME}_PATCH_LEVEL 1)
set(PACKAGE_VERSION ${VERSION})


# Revision No.
set(NPT_REVISION "20160706_1000")

set(CMAKE_CXX_COMPILER ${NPT_CXX})
message(NPT_CXX: "${NPT_CXX}")
message(CMAKE_CXX_COMPILER: "${CMAKE_CXX_COMPILER}")


#Type of REAL

if(with_real STREQUAL double)
        set(REAL_OPT "-D_REAL_IS_DOUBLE_")
endif()


# Special flags
set(NPT_LIB "Npatch")


configure_file( config.h.in.cmake config.h )
configure_file( npt-config.in.cmake npt-config)
configure_file( include/npt_Version.h.in ${PROJECT_BINARY_DIR}/include/npt_Version.h)

add_subdirectory(src)
add_subdirectory(doc)


# check
include(CheckFunctionExists)
include(CheckIncludeFiles)

check_include_files(inttypes.h HAVE_INTTYPES_H)
check_function_exists(malloc HAVE_MALLOC)
check_include_files(memory.h HAVE_MEMORY_H)
check_function_exists(sqrt HAVE_SQRT)
check_include_files(stdint.h HAVE_STDINT_H)
check_include_files(stdlib.h HAVE_STDLIB_H)
check_include_files(strings.h HAVE_STRINGS_H)
check_include_files(string.h HAVE_STRING_H)
check_function_exists(strstr HAVE_STRSTR)
check_include_files(sys/stat.h HAVE_SYS_STAT_H)
check_include_files(sys/time.h HAVE_SYS_TIME_H)
check_include_files(sys/types.h HAVE_SYS_TYPES_H)
check_include_files(unistd.h HAVE_UNISTD_H)


# install files
INSTALL(FILES ${PROJECT_BINARY_DIR}/npt-config DESTINATION ${PROJECT_NAME}/bin PERMISSIONS
OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ
WORLD_EXECUTE )
INSTALL(FILES  AUTHORS COPYING ChangeLog README LICENSE DESTINATION
 ${PROJECT_NAME}/share)
