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

rm -f "$3".so

if [ $1 = "opt" -o  $1 = "dbg" ];then
    moris_create_shared_object.sh $1 $3
    
    if [ ! -f "$3".so ];then
        exit
    fi
    moris_run.sh $1 $2 $3 $4
else
    echo ""
    echo " Error - incorrect option"
    echo ""
fi
