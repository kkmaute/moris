#!/bin/bash

### if shell script not accepted by unix, do:
#   chmod +x compile_size_scaling.sh
#   sed -i -e 's/\r$//' compile_size_scaling.sh

### to run script:
#   ./compile_size_scaling.sh >& log_compile_size_scaling.txt &
#   tail -f compile_size_scaling.txt

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
ScaleNames=("64" "32" "16" "8" "4" "2" "1")
XScales=("8" "4" "4" "2" "2" "1" "1")
YScales=("8" "8" "4" "4" "2" "2" "1")
#        1x1 2x1 2x2 4x2 4x4 8x4 8x8


for i in ${!ScaleNames[@]}; do
    
    # make sure to start from the base directory
    cd $BaseDirectory

    # get the new folder name
    NewDirName="Size_1_over_${ScaleNames[$i]}"
    
    # output
    echo ""
    echo "======================================================"
    echo "Prepare run for $NewDirName"
    
    # delete the directory if it already exists
    if [ -d  $NewDirName ]; then
        echo "Deleting old directory $NewDirName"
        rm -rf $NewDirName 
    fi
    
    # create directory for run
    echo "Creating new directory $NewDirName"
    mkdir $NewDirName
    
    # go to the new directory
    cd $NewDirName
    
    # create sub-directories for output
    for k in ${!SubDirectoriesToCreate[@]}; do
        mkdir "${SubDirectoriesToCreate[$k]}"
    done
    echo "Created sub-directories for output."
    
    # copy the input file
    cp $BaseFilePath ./
    echo "Copied template input file into new directory."
    
    # modify the input file for the current run
    echo -n "Changing parameters ... "
    ReplaceString="s/i-0000/${XScales[$i]}/1"
    sed -i $ReplaceString $FileName.cpp
    ReplaceString="s/i-0001/${YScales[$i]}/1"
    sed -i $ReplaceString $FileName.cpp
    echo "Done." 
    
    # compile the input file
    echo -n "Compiling input file ... "
    $MORISROOT/share/scripts/create_shared_object.sh . build_opt $FileName >& log_compile.txt
    echo "Done."
    
    # run it
    # echo -n "Perform simulation ... "
    # srun -n 4 --export=LD_LIBRARY_PATH,PWD $MRO ./$FileName.so >& log.txt
    # echo "Done."

done

#############################################################
echo ""
echo "======================================================"
echo "============= Done, all runs completed. =============="
echo "======================================================"
#############################################################