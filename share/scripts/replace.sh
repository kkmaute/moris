#!/bin/sh
#====================================================================================
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#====================================================================================
if [ ! "$1" ];then
  echo ""
  echo "=============================================================="
  echo ""
  echo " replace string in file(s)"
  echo ""
  echo " input: <file/file identifier> <old string> <new string>"
  echo ""
  echo "=============================================================="
  echo ""
fi

# read input data

if [ ! "$1" ];then
  echo "file or identifier"
  read fid
else
  fid=$1
fi

if [ ! "$2" ];then
  echo "old string"
  read old
else
  old=$2
fi

if [ ! "$3" ];then
  echo "new string"
  read new
else
  new=$3
fi

# loop over all files

liste=`ls "$fid"`

#liste=$fid

for file in $liste;do

echo "... processing $file"

ed -s $file <<..EOFF
    g/$old/s/$old/$new/g
    w
..EOFF

done

exit

