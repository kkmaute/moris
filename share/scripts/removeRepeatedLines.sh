#!/bin/sh

REFSRC="."

rm -f fileList

find $REFSRC -regex '.*.\.[f,c,h,C,l,y,m]' -print >>  fileList

find $REFSRC -regex '.*inc' -print >>  fileList

find $REFSRC -regex '.*hpp' -print >>  fileList

find $REFSRC -regex '.*cpp' -print >>  fileList

find $REFSRC -regex '.*tpp' -print >>  fileList

find $REFSRC -regex '.*cmake' -print >>  fileList

files=`cat fileList`

for file in $files;do

  echo ""
  echo "=============================================================="
  echo ""
  echo " processing $file"
  echo ""
  echo "=============================================================="
  echo ""

  RC=1  
  cp $file tmp1

  while [ "$RC" -gt 0 ];do
   
      cat $file | awk 'BEGIN{n=0}{ if (NF == 0) { n=n+1; if ( n == 2 ){n=n-1}else{print $0}}else{n=0;print $0}}' > tmp2

       diff tmp1 tmp2 > /dev/null

       RC=$?
 
       cp tmp2 tmp1
              
  done
  
  cp tmp2 $file

done

exit  
