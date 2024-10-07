#!/bin/sh

# file to get total number of clones and total number of unique users cloning in 2 week window
#
# adjust start and end dates below

startday=01
startmonth=05
startyear=23

endday=01
endmonth=06
endyear=24

list=`ls *-*-*` 

nclones=0
nunique=0
nsets=0

for set in $list;do

   n=`echo $set | awk -F '-' -v sm=$startmonth -v sd=$startday -v sy=$startyear -v em=$endmonth -v ed=$endday -v ey=$endyear '{ 
          n=0;
          shash=sd+12*sm+365*sy
          ehash=ed+12*em+365*ey
          chash=$1+12*$2+365*$3
          if ( chash >= shash ) {
             if ( chash <= ehash ) { n=1 }
          };
          print n }'`
                    
   if [ $n -eq 1 ];then      
            

   nsets=`expr $nsets + 1`
      
   clones=`cat $set | awk -F ':' '{print $2}' | awk -F ',' '{print $1}'`

   nclones=`expr $nclones + $clones`
   
   uniques=`cat $set | awk -F ':' '{print $3}' | awk -F ',' '{print $1}'`

   nunique=`expr $nunique + $uniques`
   
   echo "$set $nsets  $nclones  $nunique"
   
   fi
   
done   
