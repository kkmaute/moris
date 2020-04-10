#!/bin/sh
#
echo ""
echo " Script to insert line with user-defined text above a line that contains user-defined string"
echo ""
echo " Input: 1. directory name  - files in all subdirectory will be processed"
echo "        2. key string      - key to indentify line above which new line with user-defined string is inserted"
echo "        3. new line string - new string"
echo ""
echo " Example command line input: project '#include <armadillo>' '#define ARMA_ALLOW_FAKE_GCC'
echo ""
#
 
if [ $1 ];then
    dir=$1
else
    echo "provide directory name"
    exit
fi

if [ $2 ];then
    key=$2
else
    echo "no key string provided"
    exit
fi

if [ $3 ];then
    newline=$3
else
    echo "no string for line to be inserted provided"
    exit
fi
  
echo ""
echo "directory: $dir"
echo "key:       $key"
echo "new line:  $newline"

echo ""
echo " continue (y/n)"
read ans

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
