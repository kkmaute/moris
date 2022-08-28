#====================================================================================
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#====================================================================================
#!/bin/sh

XXD=kdiff3

echo ""
echo "=============================================================="
echo ""
echo " Comparison of Source Files with the following Suffixes:" 
echo ""
echo " f,c,h,C,l,y,m,inc,hpp,tex,defs,sh,cmake"
echo ""
echo " and files named:"
echo ""
echo " Makefile, makefile, runme.sh"
echo ""
echo ""
echo " Input:  directory with reference sources"
echo ""
echo "         need to edit script if $XXD not available"
echo ""
echo "=============================================================="
echo ""

chkNew='n'

if [ ! "$1" ];then
  echo "reference source directory:"
  read REFSRC
else
  REFSRC=$1
fi
  

dellist=1

if [ -f fileList ];then

  echo "use existing file list (y/n)"
  read ans

  dellist=0
  
  if [ "$ans" = 'n' ];then 
    
    dellist=1

    rm -f fileList

  fi
  
fi

if [ "$dellist" = '1' ];then

  find $REFSRC -regex '.*Makefile' -print >  fileList

  find $REFSRC -regex '.*makefile' -print >>  fileList

  find $REFSRC -regex '.*runme.sh' -print >>  fileList

  find $REFSRC -regex '.*.\.[f,c,h,C,l,y,m]' -print >>  fileList

  find $REFSRC -regex '.*inc' -print >>  fileList

  find $REFSRC -regex '.*hpp' -print >>  fileList

  find $REFSRC -regex '.*cpp' -print >>  fileList

  find $REFSRC -regex '.*tpp' -print >>  fileList

  find $REFSRC -regex '.*tex' -print >>  fileList

  find $REFSRC -regex '.*defs' -print >>  fileList

  find $REFSRC -regex '.*sh' -print >>  fileList

  find $REFSRC -regex '.*cmake' -print >>  fileList

  find $REFSRC -regex '.*txt' -print >>  fileList

  find $REFSRC -regex '.*slurm*' -print >>  fileList

  find $REFSRC -regex '.*template*' -print >>  fileList
fi


liste=`cat fileList`

 
for file in $liste;do

  bn=`echo $file | awk -F "$REFSRC/" '{print $2}'`
  kf="$bn"

  echo ""
  echo "=============================================================="
  echo ""
  echo " processing $kf"
  echo ""
  echo "=============================================================="
  echo ""


  if [ ! -f "$kf" ];then
    echo "$bn does not exist in current directory"
   if [ "$chkNew" = "y" ];then
     if [ "$ans" = "y" ];then
       echo "copy new version"
       read ans
       if [ "$ans" = "y" ];then
         mv $kf $kf".old"
         cp $file $kf
         echo "new version copied, old version saved" 
       fi
      fi
	fi
    continue
  fi

  diff --ignore-space-change $file $kf > /dev/null

  if [ "$?" != 0 ];then
   $XXD $file $kf &
   echo ""
   echo "$bn is different in new FEM"
   echo ""
   echo "do you want to change your current version (y/n)"
   read ans
   if [ "$ans" = "y" ];then
     echo ""
     echo "copy new version"
     read ans
     echo ""
     if [ "$ans" = "y" ];then
       mv $kf $kf".old"
       cp $file $kf
       echo "new version copied, old version saved" 
     else
       nedit  $kf &
       nedit -read -xrm "nedit*background: lightgreen" -xrm "nedit*foreground: red" $file &
       echo ""
       echo "press any key to continue"
       read ans
     fi
   fi
  fi

done

rm fileList

exit  

