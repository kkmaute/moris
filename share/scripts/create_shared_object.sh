#!/bin/sh
#
echo ""
echo "=============================================================="
echo ""
echo " Script to build shared object for input files"
echo ""
echo " Input: 1. name of working directory of source file and share object"
echo "        2. name of moris build directory under "'$MORISROOT'
echo "        3. name of input file (without cpp suffix)"
echo ""
echo " Output: .so file"
echo ""
echo " Example: "create_shared_object_for_input_file.sh . build myinput 
echo ""
echo "=============================================================="
echo ""

if [ ! $MORISROOT ];then
    echo " Error: MORISROOT environment variable not set"
    echo ""
    exit
fi 

if [ ! -d $MORISROOT ];then
    echo " Error: $MORISROOT does not exist"
    echo ""
    exit
fi 

if [ "$1" ];then
    workdir=$1
    
    if [ ! -d $workdir ];then
       echo " Error: $workdir does not exist"
       echo ""
       exit
    fi
else
    echo " Error: no working directory provided"
    echo ""
    exit
fi

if [ "$2" ];then
    builddir=$2

    if [ ! -d "$MORISROOT/$builddir" ];then
       echo " Error: $workdir does not exist"
       echo ""
       exit
    fi
else
    echo " Error: no moris build directory provided"
    echo ""
    exit
fi
 
if [ "$3" ];then
    cppfile=$3

    if [ ! -f "$cppfile.cpp" ];then
       echo " Error: $cppfile.cpp does not exist"
       echo ""
       exit
    fi
else
    echo " Error: no source file provided"
    echo ""
    exit
fi

cd $workdir

workdir=`pwd`

mv $MORISROOT/projects/mains/input_file.cpp /tmp/.

cp  $cppfile.cpp $MORISROOT/projects/mains/input_file.cpp

touch $MORISROOT/projects/mains/input_file.cpp

cd $MORISROOT/$builddir

make shared_object_file

cd $workdir

mv $MORISROOT/$builddir/lib/libinput_file.so $cppfile.so

rm $MORISROOT/projects/mains/input_file.cpp

mv /tmp/input_file.cpp $MORISROOT/projects/mains/.
