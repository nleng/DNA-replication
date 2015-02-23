#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.64])

# Versioning stuff
m4_define(simtools_version_major, 1)
m4_define(simtools_version_minor, 0)
m4_define(simtools_version_micro, 0)

AC_INIT([human_net], [simtools_version_major().simtools_version_minor().simtools_version_micro()], [garak@zombiepiratesfromspace.eu])

SIMTOOLS_VERSION_MAJOR=simtools_version_major()
SIMTOOLS_VERSION_MINOR=simtools_version_minor()
SIMTOOLS_VERSION_MICRO=simtools_version_micro()
AC_SUBST(SIMTOOLS_VERSION_MAJOR)
AC_SUBST(SIMTOOLS_VERSION_MINOR)
AC_SUBST(SIMTOOLS_VERSION_MICRO)

AC_CONFIG_SRCDIR([src/binary_heap.hh])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Checks for all kinds of stuff
AC_SEARCH_LIBS([strerror],[cposix])
AC_SUBST(CXXFLAGS, $CXXFLAGS)
AC_USE_SYSTEM_EXTENSIONS
AC_PROG_CXX
AC_LANG([C++])
AC_TYPE_SIZE_T
AC_DISABLE_STATIC
LT_INIT
AC_PROG_INSTALL

# Turn on warnings
[CXXFLAGS="${CXXFLAGS}"]

# Use C++11
[CXXFLAGS="${CXXFLAGS} -std=c++11"]
# Define some switches

AC_MSG_CHECKING(Wether to enable debug symbols...)
AC_ARG_ENABLE([debug], [AS_HELP_STRING([--enable-debug@<:@=full@:>@],[compile debug information into library [default=yes] ])],
    [AC_MSG_RESULT(yes)]
    [AC_DEFINE([DEBUG],1,[compile in debug info])]
    [CXXFLAGS=" ${CXXFLAGS} -ggdb3"],
    [AC_MSG_RESULT(no)]
    )	


AC_MSG_CHECKING(wether to enable optimizations...)
AC_ARG_ENABLE([optimization], [AS_HELP_STRING([--disable-optimization],[disable optimization [default=enabled] ])],
    [AC_MSG_RESULT(no)]
    ,
    [AC_MSG_RESULT(yes)]
    [CXXFLAGS="${CXXFLAGS} -O3"]
    )	
AC_MSG_CHECKING(wether to enable runtime debugging...)
AC_ARG_ENABLE([runtime-debug], [AS_HELP_STRING([--enable-runtime-debug],[enable timing code [default=yes] ])],
    [AC_MSG_RESULT(yes)]
    [CXXFLAGS="${CXXFLAGS} -DRUNTIME_DEBUG"],
    [AC_MSG_RESULT(no)]
    )

####################### Checks for libraries. #################################

# boost

BOOST_REQUIRE([1.40.0])
BOOST_MATH()

# AM_PATH_PYTHON([2.5])
# AX_PYTHON_DEVEL([>= '2.5'])
# [CPPFLAGS="${CPPFLAGS} ${PYTHON_CPPFLAGS}"]
# [LDFLAGS="${LDFLAGS} ${PYTHON_LDFLAGS} ${GSL_LDFLAGS}"]

# Python libraries
# Numpy
# AX_PYTHON_MODULE(numpy, fatal)
# NUMPY_DIR="${PYTHON_SITE_PKG}/numpy"
# AC_ARG_WITH([numpy], [AS_HELP_STRING([--with-numpy=DIR],[directory where numpy is installed
#         [default=${PYTHON_SITE_PKG}/numpy] ])],
#     NUMPY_DIR=$withval)
# [NUMPY_DIR="${NUMPY_DIR}/core/include/numpy"]
# AC_CHECK_HEADER([${NUMPY_DIR}/arrayobject.h],
#     [CPPFLAGS="${CPPFLAGS} -I${NUMPY_DIR}"],
#     [AC_MSG_ERROR([Numpy extension header not found])],
#     [#include "Python.h"])

# Scipy
# AC_PYTHON_MODULE(scipy, fatal)
# [SCIPY_DIR="${PYTHON_SITE_PKG}/scipy"]
# AC_ARG_WITH([scipy], [AS_HELP_STRING([--with-scipy=DIR],[scipy install directory
#     [default=${PYTHON_SITE_PKG}/scipy] ])],
#     [SCIPY_DIR=$withval])
# AC_CHECK_HEADER([${SCIPY_DIR}/weave/scxx/object.h],
#     [AC_DEFINE([HAVE_SCIPY], [1], [using scipy\'s weave])]
#     [CPPFLAGS="${CPPFLAGS} -I${SCIPY_DIR}"],
#     [AC_MSG_ERROR([Scipy weave/scxx header not found])])


# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
# AC_DEFINE_UNQUOTED([INSTALL_PREFIX],"${prefix}", [Python prefix])
# AC_DEFINE_UNQUOTED([PYTHON_DIR], "${PYTHON_SITE_PKG}", [The directory name for the site-packages subdirectory of the standard Python install tree.])

AC_DEFINE_UNQUOTED([CXXFLAGS],"${CXXFLAGS}", [C++ compilation options])
AC_DEFINE_UNQUOTED([CPPFLAGS],"${CPPFLAGS}", [C++ preprocessor compilation options])
AC_DEFINE_UNQUOTED([LDFLAGS],"${LDFLAGS}", [Linker options])

AC_DEFINE_UNQUOTED([AUTHOR], "Daniel Löb <garak@zombiepiratesfromspace.de>", [Program author(s)])

AC_DEFINE_UNQUOTED([COPYRIGHT], "Copyright (C) 2009-2011 Daniel Löb\nThis is free software; see the source for copying conditions. There is NO\nwarranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.", [copyright info])

AC_DEFINE([GIT_COMMIT], "esyscmd(git show | head -n 1 | sed 's/commit //' |  grep -o -e '.\{8\}' | head -n 1 |tr -d '\n')", [git HEAD commit hash])
AC_DEFINE([GIT_COMMIT_DATE], "esyscmd(git log -1 | head -n 3 | grep 'Date:' | sed s/'Date:   '// | tr -d '\n')", [git HEAD commit date])

# global CPP flags
AC_SUBST(CPPFLAGS)

# global CXX flags
AC_SUBST(CXXFLAGS)

# CPP flags for submodules
[MOD_CPPFLAGS="-I\$(srcdir) -I\$(srcdir)/.. \
-DHAVE_CONFIG_H"]
AC_SUBST(MOD_CPPFLAGS)

# default directory for submodules
# [MOD_DIR="${pythondir}/simtools"]
# AC_SUBST(MOD_DIR)

# default LIBADD flags for submodules
# [MOD_LIBADD="${PYTHON_LDFLAGS} -l${BOOST_PYTHON_LIB} ${BOOST_THREAD_LIB}"]
# AC_SUBST(MOD_LIBADD)

# Needed for typeinfo objects to work across DSO boundaries.
[MOD_LDFLAGS="-module  -avoid-version  -export-dynamic  -no-undefined\
295	-Wl,-E"]
AC_SUBST(MOD_LDFLAGS)

AC_CONFIG_FILES([
Makefile
src/Makefile
simtools.pc
])

AC_OUTPUT
