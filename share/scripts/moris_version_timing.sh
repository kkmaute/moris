#!/bin/sh
#=================================================================================
# start of user input
#=================================================================================

# overwrite MORISROOT with director of moris git version
export MORISROOT=$HOME/codes/morisGit

# build directory
here=$MORISROOT/build_opt

# file with list of examples 
exalist=$MORISROOT/share/scripts/TimingExampleList

# file with git versions to processed in additon to current one; leave empty 
# if only current git version should be checked
gitlist=$HOME/bin/CheckGitList

#=================================================================================
# end of user input
#=================================================================================

vlist=`cat $gitlist`
vlist=$vlist" "`git log | head -1 | awk '{print $2}'`

id=0

#export boxdir=/data/home/maute/work/FemdocMorisComparision/Box/moris
#export bspdir=/data/home/maute/work/FemdocMorisComparision/SphereInBox/moris

cd $here

if [ ! "$1" = "skip" ];then

    for vers in $vlist;do

        id=`expr $id + 1`
        
        tmpdate=`git log | awk -v gh=$vers  'BEGIN{n=0} 
        { 
            if( match($0,gh) > 0){n=1}
            if ( match($0,"Date") > 0 && n==1 )
            {
               split($0,a," ")
               n=0
               print a[2]" "a[3]" "a[4]" "a[5]" "a[6]
            }
        }'`
        
        date=`date -d "$tmpdate" +%m-%d-%Y`

        echo " ==============================================="

        echo " processing $id ( $vers / $date )"

        echo " ==============================================="
        
        make clean >& /dev/null

        git checkout master
        
        git pull
        
        git checkout $vers
        
        make -j 4 >& compile.$date
        
        ctest -V >& ctest.$date
        
        #cd $boxdir
        
        #moris.sh opt 1 box >& log.$date
     
        #cd bspdir
        
        #moris.sh opt 1 Bsphere >& log.$date
       
        #cd $here 
        
    done

fi  

exalist=`cat $exalist`

for exa in $exalist; do

   echo ""
   echo "================================================================"
   echo "Processing $exa"
   echo ""

   fname=$exa.timing

   rm -f $fname

   grep $exa ctest.* | grep Passed | awk -v file=$fname 'BEGIN{mid=0;mxt=10000.0;cid=0;cur=1e9}
   {
        time=$(NF-1)
        split($0,a,":") 
        split(a[1],b,".") 
        print b[2]"  "time >> file
        if ( time < mxt  ){ mid=b[2]; mxt=time}
        if ( b[2] > cid ) { cid=b[2]; cur=time}
   }
   END{print file": min-id = "mid" mxt = "mxt" cid = "cid" cur = "cur" delta = "cur-mxt" ratio [%] = "(cur-mxt)/cur*100} ' 
done

  
