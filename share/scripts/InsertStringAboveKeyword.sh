#====================================================================================
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#====================================================================================
#!/bin/sh
#
echo ""
echo "=============================================================="
echo ""
echo " Script to insert line with user-defined text above a line that contains user-defined string"
echo ""
echo " Input: 1. directory name  - files in all subdirectory will be processed"
echo "        2. key string      - key to indentify line above which new line with user-defined string is inserted"
echo "        3. new line string - new string"
echo ""
echo " Example command line input: ~/codes/moris/projects '#include <armadillo>' '#define ARMA_ALLOW_FAKE_GCC'"
echo ""
echo " Script will create a backup of each modified file; the backup is named <file>.org"
echo ""
echo " To remove all backup files use: find <directory name> | grep '\.org' | xargs rm"
echo ""
echo "=============================================================="
echo ""
#
 
if [ "$1" ];then
    dir=$1
else
    echo " Error: no directory name provided"
    echo ""
    exit
fi

if [ "$2" ];then
    key=$2
else
    echo " Error no key string provided"
    echo ""
    exit
fi

if [ "$3" ];then
    newline=$3
else
    echo " Error: no string for line to be inserted provided"
    echo ""
    exit
fi
  
echo ""
echo " directory: $dir"
echo " key:       $key"
echo " new line:  $newline"

echo ""
echo " continue (y/n)"
read ans
echo ""

if [ ! $ans = "y" ];then
   exit
fi
 
rm -f abc
rm -f dudel 
rm -f sedscript

echo "\/$key/i $newline" > sedscript
 
find $dir | grep '\.cpp' | xargs grep "$key" > dudel
find $dir | grep '\.hpp' | xargs grep "$key" >> dudel
 
list=`cat dudel`

if [ ! $list ];then
   echo " no files matching key found"
   echo ""
   exit
fi
 
rm -f dudel
 
for obj in $list;do
 
    fname=`echo $obj | awk -F ":" '{ if ( NF > 1) {print $1}}'`
    
    if [ $fname ];then
       echo " ... processing $fname" 
        
       sed -f sedscript $fname > abc
              
       cp $fname $fname.org
       
       mv abc $fname
 
    fi
done
 
rm -f abc
rm -f dudel 
rm -f sedscript

