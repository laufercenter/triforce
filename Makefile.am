AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

lib_LTLIBRARIES = libtriforce-@TRIFORCE_API_VERSION@.la

libtriforce_@TRIFORCE_API_VERSION@_la_SOURCES = triforce/tessellation.cc \
						triforce/molecule.cc \
						triforce/spline.cc \
						triforce/datafiledigest.cc

libtriforce_@TRIFORCE_API_VERSION@_la_LDFLAGS = -version-info $(TRIFORCE_SO_VERSION)

triforce_includedir = $(includedir)/triforce-$(TRIFORCE_API_VERSION)
nobase_triforce_include_HEADERS = triforce.h \
					triforce/tessellation.h \
					triforce/molecule.h \
					triforce/spline.h \
					triforce/datafiledigest.h

triforce_libincludedir = $(libdir)/triforce-$(TRIFORCE_API_VERSION)/include
nodist_triforce_libinclude_HEADERS = triforceconfig.h

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = triforce-$(TRIFORCE_API_VERSION).pc

dist_noinst_SCRIPTS = autogen.sh

dist_doc_DATA = dat
                                            
