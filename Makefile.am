AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

lib_LTLIBRARIES = libtriforce-@TRIFORCE_API_VERSION@.la

libtriforce_@TRIFORCE_API_VERSION@_la_SOURCES = triforce/tessellation.cc \
						triforce/molecule.cc \
						triforce/data3d.cc \
						triforce/datafile.cc \
						triforce/topology.cc \
						triforce/interpolation.cc \
						triforce/integrator.cc

libtriforce_@TRIFORCE_API_VERSION@_la_LDFLAGS = -version-info $(TRIFORCE_SO_VERSION)

triforce_includedir = $(includedir)/triforce-$(TRIFORCE_API_VERSION)
nobase_triforce_include_HEADERS = triforce.h \
					triforce/tessellation.h \
					triforce/molecule.h \
					triforce/data3d.h \
					triforce/datafile.h \
					triforce/topology.h \
					triforce/interpolation.h \
					triforce/integrator.h

triforce_libincludedir = $(libdir)/triforce-$(TRIFORCE_API_VERSION)/include
nodist_triforce_libinclude_HEADERS = triforceconfig.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = triforce-$(TRIFORCE_API_VERSION).pc

dist_noinst_SCRIPTS = autogen.sh

                                            
