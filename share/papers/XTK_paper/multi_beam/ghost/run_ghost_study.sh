#!/bin/bash

### if shell script not accepted by unix, do:
#   chmod +x run_ghost_study.sh
#   sed -i -e 's/\r$//' run_ghost_study.sh

### to run script:
#   ./run_ghost_study.sh >& log_ghost_study.txt &
#   tail -f log_ghost_study.txt

#############################################################
echo ""
echo "======================================================"
echo "======== Write, Compile, and Run all files... ========"
echo "======================================================"
echo ""
#############################################################

BaseDirectory=$PWD
echo "Base directory = $BaseDirectory "

export FileName="Multi_Beam"

BaseFilePath="${BaseDirectory}/${FileName}.cpp"

### Parameters to iterate over:
GhostTitle="Ghost"
GhostNames=("On" "Off")
GhostValues=("true" "false")
# GhostNames=("On")
# GhostValues=("true")

### OffSetTitle="OffSet"
### OffSetNames=("000" "001" "002" "005" "010" "020" "030" "040" "050")
### OffSetValues=("0.0" "0.01" "0.02" "0.05" "0.1" "0.2" "0.3" "0.4" "0.5")

OffSetTitle="OffSet"
OffSetNames=("040" "050" "060" "061" "062" "065" "070" "080" "090" "095" "098" "099" "100" "110" "120")
OffSetValues=("0.4" "0.5" "0.6" "0.61" "0.62" "0.65" "0.7" "0.8" "0.9" "0.95" "0.98" "0.99" "1.0" "1.1" "1.2")
# OffSetNames=("-140" "-139" "-120" "-101" "-100" "-080" "-060" "-040" "-039" "-020" "-001" "000" "020" "040" "060" "061" "080" "099" "100" "120" "140" "160" "161" "180" "199" "200")
# OffSetValues=("-1.4" "-1.39" "-1.2" "-1.01" "-1.0" "-0.8" "-0.6" "-0.4" "-0.39" "-0.2" "-0.01" "0.0" "0.2" "0.4" "0.6" "0.61" "0.8" "0.99" "1.0" "1.2" "1.4" "1.6" "1.61" "1.8" "1.99" "2.0")

# 
for i in ${!GhostValues[@]}; do
    
    # make sure to start from the base directory
    cd $BaseDirectory

    # get the new folder name
    NewDirName="${GhostTitle}_${GhostNames[$i]}"
    
    # output
    echo ""
    echo "======================================================"
    echo "Prepare runs for $NewDirName"
    
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
    
    # process sub-directories 
    for j in ${!OffSetValues[@]}; do
    
        # get the new sub-folder name
        NewSubDirName="${OffSetTitle}_${OffSetNames[$j]}"
        
        # output
        echo ""
        echo "----------------------------------------------"
        
        # delete the sub-directory if it already exists
        if [ -d  $NewSubDirName ]; then
            echo "Deleting old sub-directory $NewSubDirName"
            rm -rf $NewSubDirName 
        fi
        
        # create sub-directory for run
        echo "Creating new sub-directory $NewSubDirName"
        mkdir $NewSubDirName
        
        # go to the new sub-directory
        cd $NewSubDirName
        
        # copy the template input file
        cp $BaseFilePath ./
        echo "Copied template input file into new directory."
        
        # modify the input file for the current run
        echo -n "Changing parameters ... "
        ReplaceString="s/i-0000/${GhostValues[$i]}/1"
        sed -i $ReplaceString $FileName.cpp
        ReplaceString="s/i-0001/${OffSetValues[$j]}/1"
        sed -i $ReplaceString $FileName.cpp
        echo "Done." 
        
        # compile the input file
        echo -n "Compiling input file ... "
        $MORISROOT/share/scripts/create_shared_object.sh . build_opt $FileName >& log_compile.txt
        echo "Done."
        
        # run it
        echo -n "Perform simulation ... "
        $MRO ./$FileName.so >& log.txt
        echo "Done."
        
        # go one back down in directory structure before processing next sub-directory
        cd ..
    
    # end: loop over offset values
    done
    
# end: loop over ghost values
done

#############################################################
echo ""
echo "======================================================"
echo "============= Done, all runs completed. =============="
echo "======================================================"
#############################################################
