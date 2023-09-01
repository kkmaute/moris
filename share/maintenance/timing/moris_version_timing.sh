#!/bin/sh
#=================================================================================
#
# script to perform timing test for moris 
#
# 1. copy this file to your $HOME
# 2. edit information in section below
# 3. run script with no other major job or process running on your machine
# 
# to just process existing timing results use: moris_version_timing.sh "skip"
#
# start of user input
#=================================================================================

# overwrite MORISROOT with director of moris git version
export MORISROOT=$HOME/codes/moris

# build directory
here=$MORISROOT/build_opt

# branch: master or xtk_refactor or main
branch=main

# file with list of examples 
exalist=$MORISROOT/share/maintenance/timing/TimingExampleList

# file with git versions to processed in additon to current one; leave empty 
# if only current git version should be checked
gitlist=$MORISROOT/share/maintenance/timing/CheckGitList_github
#gitlist=

#=================================================================================
# end of user input
#=================================================================================

cd $MORISROOT

git checkout $branch

if [ $gitlist ];then
    vlist=`cat $gitlist`
    vlist=$vlist" "`git log | head -1 | awk '{print $2}'`
else
    vlist=`git log | head -1 | awk '{print $2}'`
fi

id=0

cd $here

if [ ! -d "TimingResults" ];then
    mkdir TimingResults
fi    

if [ ! "$1" = "skip" ];then

    for vers in $vlist;do

        id=`expr $id + 1`

        make clean >& /dev/null

        git checkout $branch
        
        git pull
        
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
                
        git checkout $vers
        
        rm -r -f cmake CMakeCache.txt CMakeDoxyfile.in CMakeDoxygenDefaults.cmake
        rm -r -f CMakeFiles cmake_install.cmake CTestTestfile.cmake generated
        rm -r -f lib Makefile share

        cmake -DBUILD_ALL=ON -DMORIS_USE_EXAMPLES=ON ..  >& /dev/null
        cmake -DBUILD_ALL=ON -DMORIS_USE_EXAMPLES=ON ..  >& /dev/null
        
        make -j 4 >& TimingResults/compile.$date
        
        ctest -V >& TimingResults/ctest.$date
        
    done

fi  

cd TimingResults

exalist=`cat $exalist`

ctestlist=`ls -t -r  ctest.*`

for exa in $exalist; do

   echo ""
   echo "================================================================"
   echo "Processing $exa"
   echo ""

   fname=$exa.timing

   rm -f $fname

   # important: keep whitespace before $exa
   grep " $exa" $ctestlist | grep Passed | awk -v file=$fname 'BEGIN{mxt=10000.0;curnum=0}
   {
        time=$(NF-1)
        split($0,a,":") 
        split(a[1],b,".") 
        split(b[2],c,"-")
        numdate=10000*(c[3]-2000)+100*c[1]+c[2]
        print numdate"  "time >> file
        if ( time    < mxt    ){ mid=b[2]; mxt=time}
        if ( numdate > curnum ){ cid=b[2]; cur=time; curnum=numdate}
   }
   END{print file": min-id = "mid" mxt = "mxt" cid = "cid" cur = "cur" delta = "cur-mxt" ratio [%] = "(cur-mxt)/(1e-12+cur)*100} ' 
done

  
