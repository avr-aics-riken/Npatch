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

lib_LIBRARIES = libNpatch.a

libNpatch_a_CXXFLAGS = @REAL_OPT@ -I../include -I$(top_srcdir)/include

libNpatch_a_SOURCES = Npt.cxx FNpt.cxx

# npt_Version.h.in -> BUILD_DIR/include/npt_Version.h -> install
NPT_includedir = $(includedir)
NPT_include_HEADERS = \
   ../include/npt_Version.h

# SRC_DIR/include/* -> install
nobase_include_HEADERS = \
   ../include/Npt.h \
   ../include/FNpt.h \
   ../include/CalcGeo.h \
   ../include/CalcGeo_Matrix.h

