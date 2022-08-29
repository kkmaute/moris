#!/bin/sh
#====================================================================================
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#====================================================================================
#
# script to compile and run moris in either debug or optimized mode
#
# requires environment variabes MORISBUILDDBG and MORISBUILDOPT 
#
# input: 1 - "dbg" or "opt"
#        2 - number of processors: 1,2,3,4....
#        3 - basename of input file (without .cpp)
#        4 - run prefix, e.g., valgrind

sofile="$3.so"
ofile="$3.o"

rm -f $sofile $ofile

if [ $1 = "opt" -o  $1 = "dbg" ];then
    moris_create_shared_object.sh $1 $3
    
    if [ ! -f $sofile ];then
        echo ""
        echo " compilation failed"
        echo ""
        exit
    fi
    moris_run.sh $1 $2 $3 $4
else
    echo ""
    echo " Error - incorrect option"
    echo ""
fi

