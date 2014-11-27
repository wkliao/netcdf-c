## This is an automake file, part of Unidata's netCDF package.
# Copyright 2005, see the COPYRIGHT file for more information.

# This automake file is in charge of building the libsrc directory,
# which contains the classic library code.

include $(top_srcdir)/lib_flags.am

libnetcdf5_la_CPPFLAGS = ${AM_CPPFLAGS} -I$(top_srcdir)/libsrc

# This is the code for a dispatch table for pnetcdf
# (which has CDF5 as its magic number)

# Turn on a pre-processor flag when building a DLL for windows.
if BUILD_DLL
libnetcdf5_la_CPPFLAGS += -DDLL_EXPORT
endif # BUILD_DLL

# These files comprise the pnetcdf dispatch library code.
libnetcdf5_la_SOURCES = nc5dispatch.c \
nc5v1hpg.c nc5putget.c nc5attr.c nc5dispatch.c nc5internal.c \
nc5var.c nc5dim.c nc5ncx.c nc5ncx.h

noinst_LTLIBRARIES = libnetcdf5.la

man_MANS =

# These files are cleaned on developer workstations (and then rebuilt
# with m4), but they are included in the distribution so that the user
# does not have to have m4.
MAINTAINERCLEANFILES = $(man_MANS) nc5attrx.c nc5putgetx.c
EXTRA_DIST = nc5attr.m4 nc5ncx.m4 nc5putget.m4 $(man_MANS) CMakeLists.txt

# This tells make how to turn .m4 files into .c files.
.m4.c:
	m4 $(AM_M4FLAGS) $(M4FLAGS) $< >$@

