#!/bin/sh
 
if [ $1 ];then
    dir=$1
else
    echo "provide directory name"
    exit
fi
 
key='#include <armadillo>'
newline='#define ARMA_ALLOW_FAKE_GCC'
 
echo ""
echo "key:      $key"
echo "new line: $newline"

echo ""
echo " continue (y/n)"
read ans

if [ ! $ans = "y" ];then
   exit
fi
 
echo "\/$key/i $newline" > sedscript
 
rm -f dudel
 
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
