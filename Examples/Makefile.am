###################################################################################
#
# NpatchLib - Nagata Patch Library
#
# Copyright (c) 2015 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
###################################################################################


noinst_PROGRAMS = stl2npt


DISTCLEANFILES=*~ stl2npt
CLEANFILES=

EXTRA_DIST= README


stl2npt_SOURCES  =  stl2npt.cxx
stl2npt_CXXFLAGS = -I$(top_srcdir)/include

stl2npt_LDADD = \
-L$(top_builddir)/src -lNpatch


dist_noinst_DATA=

.PHONY:$(dist_noinst_DATA)
