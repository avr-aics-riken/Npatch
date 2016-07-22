#! /bin/sh

prefix=@prefix@
exec_prefix=@exec_prefix@
includedir=@includedir@
libdir=@libdir@

usage()
{
cat <<EOF

Usage: npatch-config [OPTION]

Known values for OPTION are:

--cxx       print C++ compiler command
--cflags    print C/C++ pre-processor and compiler flags
--libs      print library linking information for C++ program
--help      display this help and exit
--version   output version information

EOF

exit $1
}

if test $# -eq 0; then
usage 1
fi

cflags=false
libs=false

while test $# -gt 0; do
case "$1" in
-*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
*) optarg= ;;
esac

case "$1" in
--version)
cat <<EOF

Npatch - Nagata Patch Library  Version : @VERSION@ : @NPT_REVISION@

Copyright (c) 2015-2016 Advanced Institute for Computational Science, RIKEN.
All rights reserved.

EOF
exit 0
;;

--help)
usage 0
;;

--cxx)
echo @NPT_CXX@
;;

--cflags)
echo -I@NPT_DIR@/include
;;

--libs)
echo -L@NPT_DIR@/lib -l@NPT_LIB@
;;

*)
usage
exit 1
;;
esac
shift
done

exit 0
