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

SUBDIRS = @MAKE_SUB_DIRS@

bin_SCRIPTS=npt-config
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = npt-config.in LICENSE BUILD_DIR/.gitkeep \
             doc/doxygen/Doxyfile Examples
dist_data_DATA = AUTHORS COPYING ChangeLog README LICENSE NEWS README.md

clean-local:
	find ./ -name *.obj | xargs rm -f
	find ./ -name *.exe | xargs rm -f
	find ./ -name *.lib | xargs rm -f

distclean-local: clean-local
	rm -f ./bin/*
	rm -f ./lib/*
