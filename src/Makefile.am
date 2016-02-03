##############################################################################
#
# Npatchlib - Nagata Patch Library
#
# Copyright (c) 2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
##############################################################################

AUTOMAKE_OPTIONS = subdir-objects

lib_LIBRARIES = libNpatch.a
libNpatch_a_CXXFLAGS = -I$(top_srcdir)/include
libNpatch_a_SOURCES = \
     Npt.cxx \
     FNpt.cxx

nobase_include_HEADERS = \
  ../include/CalcGeo.h \
  ../include/CalcGeo_Matrix.h \
  ../include/FNpt.h \
  ../include/Npt.h


EXTRA_DIST = Makefile_hand depend.inc
