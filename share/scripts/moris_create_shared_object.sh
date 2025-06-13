#!/bin/sh

#====================================================================================
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#====================================================================================

#
# script to compile moris input file in either debug or optimized mode
#
# requires environment variabes MORISBUILDDBG and MORISBUILDOPT 
#
# input: 1 - "dbg" or "opt"
#        2 - basename of input file (without .cpp)
#        3 - keyword dbg to keep link to input file

if [ ! "$MORISBUILDDBG" -o ! "$MORISBUILDOPT" ];then
   echo ""
   echo " Error - environment variabes MORISBUILDDBG and MORISBUILDOPT required"
   echo ""
   exit
fi

if [ $1 = "opt" ];then
    echo ""
    echo " ... compiling against optimized moris"
    echo ""
    $MORISROOT/share/scripts/create_shared_object.sh . $MORISBUILDOPT $2 $3
else
   if [ $1 = "dbg" ];then
       echo ""
       echo " ... compiling against debug moris"
       echo ""
       $MORISROOT/share/scripts/create_shared_object.sh . $MORISBUILDDBG $2 $3
   else
       echo ""
       echo " Error - incorrect option. Compile with 'opt' or 'dbg'. as the first argument"
       echo ""
   fi
fi

