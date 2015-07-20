###
###################################################################################
#
# NpatchLib - Nagata Patch Library
#
# Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
###################################################################################
###

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])

AC_INIT([NpatchLib],[0.8.1],[keno@riken.jp],[NpatchLib])

echo "----------------------------------------"
echo "Running NpatchLib Configure Script"
echo "----------------------------------------"


AC_CONFIG_SRCDIR([src/npt_Data.cxx])


#
# Specify directory for auxillary build tools (e.g., install-sh, config.sub, config.guess) and M4 files.
#

AC_CONFIG_AUX_DIR(config)


#
# To prevent the source codes from breaking, we suppress to do configure in the source directory.
#

if test -e configure.ac ; then
echo "Please run configure from a separate build directory. It is not allowd to run configure in the source directory."
echo "Pleaes try ..."
echo " "
echo "$ cd BUILD_DIR"
echo "$ ../configure [options]"
exit
fi


#
# Configure should create config.h from config.h.in
#

AC_CONFIG_HEADERS([config.h])


#
# Use automake
#

AM_INIT_AUTOMAKE()


#
# revision No.
#

NPT_REVISION="20150720_1800"
AC_SUBST(NPT_REVISION)


#
# default prefix setting
#

AC_PREFIX_DEFAULT([/usr/local/NpatchLib])
if test "$prefix" = "NONE" ; then
NPT_INST_DIR=/usr/local/NpatchLib
else
NPT_INST_DIR="$prefix"
fi


NPT_CXX="$CXX"
AC_SUBST(NPT_CXX)


#
# example
#
AC_ARG_WITH(example, [AC_HELP_STRING([--with-example=(no|yes)],[make example [no]])], , with_example=no)



#
# Npatch special flags
#

NPT_CFLAGS="-I$NPT_INST_DIR/include"
AC_SUBST(NPT_CFLAGS)

NPT_LDFLAGS="-L$NPT_INST_DIR/lib"
AC_SUBST(NPT_LDFLAGS)

NPT_LIBS="-lNpatch"
AC_SUBST(NPT_LIBS)



#
# MAKE_SUB_DIRS
#

AC_SUBST(MAKE_SUB_DIRS)
MAKE_SUB_DIRS="src doc"

if test x"$with_example" = x"yes" ; then
MAKE_SUB_DIRS="$MAKE_SUB_DIRS"" Examples"
fi


#
# Checks for programs.
#

AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB



#
# Checks for header files.
#

AC_CHECK_HEADERS([stdlib.h string.h sys/timeb.h])


#
# Checks for typedefs, structures, and compiler characteristics.
#

AC_CHECK_HEADER_STDBOOL
AC_C_INLINE


#
# Checks for library functions.
#

AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt strstr])



AC_CONFIG_FILES([Makefile \
                 src/Makefile \
                 Examples/Makefile \
                 doc/Makefile \
                 npt-config \
                 include/npt_Version.h \
                 ])

AC_OUTPUT


chmod +x ./npt-config

echo "---------------------------------------------"
echo "Finished Running NpatchLib Configure Script"
echo "---------------------------------------------"
