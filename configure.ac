###
###################################################################################
#
# NpatchLib - Nagata Patch Library
#
# Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
###################################################################################
###

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.57])

AC_INIT([NpatchLib],[1.0.0],[keno@riken.jp],[NpatchLib])

echo "----------------------------------------"
echo "Running NpatchLib Configure Script"
echo "----------------------------------------"


AC_CONFIG_SRCDIR([src/Npt.cxx])


#
# Specify directory for auxillary build tools (e.g., install-sh, config.sub, config.guess) and M4 files.
#

AC_CONFIG_AUX_DIR(config)


#
# To prevent the source codes from breaking, we suppress to do configure in the source directory.
#

if test -e configure.ac ; then
echo "Please run configure from a separate build directory. It is not allowed to run configure in the source directory."
echo "Please try ..."
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

NPT_REVISION="20160215_1000"
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
# Npatch special flags
#

NPT_CFLAGS="-I$NPT_INST_DIR/include"
AC_SUBST(NPT_CFLAGS)

NPT_LDFLAGS="-L$NPT_INST_DIR/lib"
AC_SUBST(NPT_LDFLAGS)

NPT_LIBS="-lNpatch"
AC_SUBST(NPT_LIBS)

#
# Type of REAL
#
AC_ARG_WITH(real, [AC_HELP_STRING([--with-real=(float|double)],[Specify real type [float]])], , with_real=float)
AC_SUBST(REAL_OPT)

if test x"$with_real" = x"double" ; then
  REAL_OPT=-D_REAL_IS_DOUBLE_
  NPT_CFLAGS="$NPT_CFLAGS"" $REAL_OPT"
else
  REAL_OPT=
fi


#
# MAKE_SUB_DIRS
#

AC_SUBST(MAKE_SUB_DIRS)
MAKE_SUB_DIRS="src doc"


#
# Checks for programs.
#

AC_PROG_CXX
AC_PROG_RANLIB



#
# Checks for header files.
#

AC_CHECK_HEADERS([stdlib.h string.h sys/time.h])


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
                 doc/Makefile \
                 npt-config \
                 include/npt_Version.h \
                 ])

AC_OUTPUT


chmod +x ./npt-config

echo "---------------------------------------------"
echo "Finished Running NpatchLib Configure Script"
echo "---------------------------------------------"
