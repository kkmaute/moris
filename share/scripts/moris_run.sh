#!/bin/sh

#====================================================================================
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#====================================================================================

#
# script to run moris in either debug or optimized mode
#
# requires environment variabes MORISBUILDDBG and MORISBUILDOPT 
#
# input: 1 - "dbg" or "opt"
#        2 - number of processors: 1,2,3,4....
#        3 - shared object file (without .so)
#        4 - run prefix, e.g., valgrind

if [ ! "$MORISBUILDDBG" -o ! "$MORISBUILDOPT" ];then
   echo ""
   echo " Error - environment variabes MORISBUILDDBG and MORISBUILDOPT required"
   echo ""
   exit
fi

if [ $1 = "opt" ];then
    echo ""
    echo " ... running optimized moris on $2 processors"
    echo ""
    mpirun -np $2 $4 $MORISROOT/$MORISBUILDOPT/projects/mains/moris "$3.so"
else
   if [ $1 = "dbg" ];then
       echo ""
       echo " ... running debug moris on $2 processors"
       echo ""
       mpirun -np $2 $4 $MORISROOT/$MORISBUILDDBG/projects/mains/moris "$3.so"
   else
       echo ""
       echo " Error - incorrect option"
       echo ""
   fi
fi

