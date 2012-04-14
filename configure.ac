## TRIFORCE

AC_INIT([TRIFORCE], [0.39], [nils.drechsel@upf.edu], [TRIFORCE], [http://cbbl.imim.es])
AC_PREREQ([2.59])
AM_INIT_AUTOMAKE([1.10 -Wall no-define])

AC_CONFIG_HEADERS([config.h triforceconfig.h])
AC_PROG_CXX
LT_INIT([disable-static])

AC_SUBST([TRIFORCE_SO_VERSION], [1:0:0])
AC_SUBST([TRIFORCE_API_VERSION], [0.39])

AC_CONFIG_FILES([Makefile
		                  triforce-${TRIFORCE_API_VERSION}.pc:triforce.pc.in])
AC_OUTPUT
~      
