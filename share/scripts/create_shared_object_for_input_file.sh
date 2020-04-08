#!/bin/sh
 
dir=$1
 
mrs_dir=$2

build=$3

abc=$4

mv $mrs_dir/projects/mains/input_file.cpp /tmp/.

ln -s $abc.cpp $mrs_dir/projects/mains/input_file.cpp

cd $build

make input_file

cd $dir

cp $build/lib/libinput_file.so $abc.so

rm $mrs_dir/projects/mains/input_file.cpp

mv /tmp/input_file.cpp $MORISROOT/projects/mains/.