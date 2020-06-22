#!/bin/sh
#
# script to compile and run moris in either debug or optimized mode
#
# requires environment variabes MORISBUILDDBG and MORISBUILDOPT 
#
# input: 1 - "dbg" or "opt"
#        2 - number of processors: 1,2,3,4....
#        3 - basename of input file (without .cpp)
#        4 - run prefix, e.g., valgrind

if [ $1 = "opt" -o  $1 = "dbg" ];then
    moris_create_shared_object.sh $1 $3
    moris_run.sh $1 $2 $3 $4
else
    echo ""
    echo " Error - incorrect option"
    echo ""
fi
