#!/bin/sh
# After running this script, the package will be ready to configure/make/make install.

top_srcdir=`dirname $0`
test -z "$top_srcdir" && top_srcdir=.
top_builddir=$top_srcdir
echo $top_srcdir

PROJECT="replication"

(test -f $top_srcdir/configure.ac) || {
    echo -n "ERROR: Directory \"$srcdir\" does not seem to be the top-level"
    echo -n " package directory."
    exit 1
}

aclocal -I m4 || exit $?
autoheader || exit $?
libtoolize -f || exit $?
automake --add-missing --copy || exit $?
autoconf || exit $?