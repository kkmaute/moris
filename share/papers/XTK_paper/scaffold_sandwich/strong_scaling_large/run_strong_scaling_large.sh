#!/bin/bash

### if shell script not accepted by unix, do:
#   chmod +x run_strong_scaling.sh
#   sed -i -e 's/\r$//' run_strong_scaling.sh

### to run script:
#   ./run_strong_scaling.sh >& log_strong_scaling.txt &
#   tail -f run_strong_scaling.txt

#############################################################
echo ""
echo "======================================================"
echo "======== Write, Compile, and Run all files... ========"
echo "======================================================"
echo ""
#############################################################

BaseDirectory=$PWD
echo "Base directory = $BaseDirectory "

export FileName="nTop_Sandwich"

BaseFilePath="${BaseDirectory}/${FileName}.cpp"

SubDirectoriesToCreate=("XTK")

### For Alpine:
NumProcs=("1" "2" "4" "8" "16" "32" "64")
ProcSplits=("1,1,1" "2,1,1" "2,2,1" "4,2,1" "4,4,1" "8,4,1" "8,8,1")


for i in ${!NumProcs[@]}; do
    
    # make sure to start from the base directory
    cd $BaseDirectory

    # get the new folder name
    NewDirName="np_${NumProcs[$i]}"
    
    # output
    echo ""
    echo "======================================================"
    echo "Run $NewDirName"
    
    # # delete the directory if it already exists
    # if [ -d  $NewDirName ]; then
    #     echo "Deleting old directory $NewDirName"
    #     rm -rf $NewDirName 
    # fi
    # 
    # # create directory for run
    # echo "Creating new directory $NewDirName"
    # mkdir $NewDirName
    
    # go to the new directory
    cd $NewDirName
    
    # # create sub-directories for output
    # for k in ${!SubDirectoriesToCreate[@]}; do
    #     mkdir "${SubDirectoriesToCreate[$k]}"
    # done
    # echo "Created sub-directories for output."
    # 
    # # copy the input file
    # cp $BaseFilePath ./
    # echo "Copied template input file into new directory."
    # 
    # # modify the input file for the current run
    # echo -n "Changing parameters ... "
    # ReplaceString="s/i-0000/${ProcSplits[$i]}/1"
    # sed -i $ReplaceString $FileName.cpp
    # echo "Done." 
    # 
    # # compile the input file
    # echo -n "Compiling input file ... "
    # $MORISROOT/share/scripts/create_shared_object.sh . build_opt $FileName >& log_compile.txt
    # echo "Done."
    
    # run it
    echo -n "Perform simulation ... "
    srun -n ${NumProcs[$i]} --export=LD_LIBRARY_PATH,PWD $MRO ./$FileName.so >& log.txt
    echo "Done."

done

#############################################################
echo ""
echo "======================================================"
echo "============= Done, all runs completed. =============="
echo "======================================================"
#############################################################
