#!/bin/bash
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

# overwrite MORISROOT with directory of moris git version
export MORISROOTORG=$HOME/codes/moris

# build directory
here=$MORISROOTORG/build_opt

# branch: master or xtk_refactor or main
branch=main

# file with list of examples 
exalist=$MORISROOTORG/share/maintenance/timing/TimingExampleList

# set current version (t13 or t15)
ctvers=t15

# resource file for current and old versions
bashrc_t13=$HOME/BASHRC_MORIS_T13
bashrc_t15=$HOME/.bashrc_moris

# file with git versions to processed in additon to current one; leave empty 
# if only current git version should be checked
gitlist=$MORISROOTORG/share/maintenance/timing/CheckGitList_github
#gitlist=

#=================================================================================
# end of user input
#=================================================================================

cd $MORISROOTORG

git checkout $branch

if [ $gitlist ];then
    vlist=`cat $gitlist | awk '{print $1}'`
    ulist=`cat $gitlist | awk '{print $2}'`
    vlist=$vlist" "`git log | head -1 | awk '{print $2}'`
    ulist=$ulist" "`echo $ctvers`
else
    vlist=`git log | head -1 | awk '{print $2}'`
    ulist=`echo $ctvers`
fi

id=0

cd $here

cp $MORISROOTORG/share/cmake/find_modules/FindOPENBLAS.cmake /tmp/FindOPENBLAS.cmake 

if [ ! -d "TimingResults" ];then
    mkdir TimingResults
fi    

if [ ! "$1" = "skip" ];then

    for vers in $vlist;do

        id=`expr $id + 1`
        
        urs=`echo $ulist | awk -v id=$id '{split($0, list);print list[id]}'`

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
                        
        doskip=0
        
        if [ "$urs" = "t13" ]; then    
            if [ ! -f $bashrc_t13 ]; then
                echo "bashrc_t13: $bashrc_t13 does not exist"
                echo "skipping timing test for $vers"
                doskip=1
            else
                echo "trilinos version 13"
                source $bashrc_t13
                echo $Trilinos_DIR
            fi
        fi
         
        if [ "$urs" = "t15" ]; then    
            if [ ! -f $bashrc_t15 ]; then
                echo "bashrc_t15: $bashrc_t15 does not exist"
                echo "skipping timing test for $vers"
                doskip=1
            else
                echo "trilinos version 15"
                source $bashrc_t15
                echo $Trilinos_DIR
            fi
        fi
        
        if [ "$doskip" = "0" ];then
        
            git checkout $vers
        
            rm -r -f cmake CMakeCache.txt CMakeDoxyfile.in CMakeDoxygenDefaults.cmake
            rm -r -f CMakeFiles cmake_install.cmake CTestTestfile.cmake generated
            rm -r -f lib Makefile share
            
            mv $MORISROOTORG/share/cmake/find_modules/FindOPENBLAS.cmake /tmp/FindOPENBLAS.cmake.current

            cp /tmp/FindOPENBLAS.cmake $MORISROOTORG/share/cmake/find_modules/FindOPENBLAS.cmake
        
            echo "MORIS cmake for $vers at $date"            > TimingResults/cmake.$date
            echo " "                                        >> TimingResults/cmake.$date
            cmake -DBUILD_ALL=ON -DMORIS_USE_EXAMPLES=ON .. >> TimingResults/cmake.$date 2>&1
            cmake -DBUILD_ALL=ON -DMORIS_USE_EXAMPLES=ON .. >> TimingResults/cmake.$date 2>&1 
        
            echo "MORIS compilation log for $vers at $date"  > TimingResults/compile.$date
            echo " "                                        >> TimingResults/compile.$date
            make -j 5                                       >> TimingResults/compile.$date 2>&1
              
            echo "MORIS ctest for $vers at $date"  > TimingResults/ctest.$date
            echo " "                              >> TimingResults/ctest.$date
            ctest -V                              >> TimingResults/ctest.$date 2>&1            
            
            mv /tmp/FindOPENBLAS.cmake.current $MORISROOTORG/share/cmake/find_modules/FindOPENBLAS.cmake
            
        fi        
    done
    
    git checkout $branch
        
    git pull
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

  
#   to extract timing data of a specific example use for example for test #86 and ctest.11-30-2022
#
#    grep '86: |  |_ElapsedTime' ctest.11-30-2022 | awk -F '=' '{print $2}' | awk '{print $1}' > 86.11-30-2022
#
