#!/bin/sh
# Exmaple command line arg:
# 		~/Codes/moris/share/scripts/create_shared_object_for_input_file.sh . /home/doble/Codes/moris /home/doble/Codes/moris/build/ /home/doble/Downloads/benchmark04
#

 
dir=$1
 
mrs_dir=$2

build=$3

abc=$4

mv $mrs_dir/projects/mains/input_file.cpp /tmp/.


ln -s $abc.cpp $mrs_dir/projects/mains/input_file.cpp

touch $mrs_dir/projects/mains/input_file.cpp

cd $build

make input_file

cd $dir

cp $build/lib/libinput_file.so $abc.so

rm $mrs_dir/projects/mains/input_file.cpp

mv /tmp/input_file.cpp $MORISROOT/projects/mains/.