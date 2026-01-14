#!/bin/bash

### if shell script not accepted by unix, do:
#   chmod +x run_convergence_study.sh
#   sed -i -e 's/\r$//' run_convergence_study.sh

### to run script:
#   ./run_convergence_study.sh >& log_convergence_study.txt &
#   tail -f log_convergence_study.txt

#############################################################
echo ""
echo "======================================================"
echo "======== Write, Compile, and Run all files... ========"
echo "======================================================"
echo ""
#############################################################

echo "PID: $$"

BaseDirectory=$PWD
echo "Base directory = $BaseDirectory "

export FileName="Multi_Beam"
export ReferenceFileName="${FileName}_Reference"
export CoarseFileName="${FileName}_Coarse"
export CompareFileName="${FileName}_Compare"

BaseFilePathReference="${BaseDirectory}/${ReferenceFileName}.cpp"
BaseFilePathCoarse="${BaseDirectory}/${CoarseFileName}.cpp"
BaseFilePathCompare="${BaseDirectory}/${CompareFileName}.cpp"

RefSolnDirName="Reference_Solution" 
CoarseSolnDirName="Coarse_Solution" 
CompareSolnDirName="Compare_Solution" 

### Parameters to iterate over:
ReferenceRefinement="7"
ReferenceElementSize="0.0078125"
EnrichmentValues=("true" "false")
RefinementValues=("0" "1" "2" "3" "4" "5" "6")
ElementSizes=("1.0" "0.5" "0.25" "0.125" "0.0625" "0.03125" "0.015625")

### Debug or optimized version of MORIS
UseDebug=0

# loop over Enrichment Values
for e in ${!EnrichmentValues[@]}; do

    # make sure to start this process from the base directory
    cd $BaseDirectory

    # get the new folder name
    EnrDirName="Enrichment_${EnrichmentValues[$e]}"
    echo ""
    echo "======================================================"
    echo "Prepare runs for $EnrDirName"
    if [ -d  $EnrDirName ]; then
        echo "Deleting old directory $EnrDirName"
        rm -rf $EnrDirName 
    fi
    echo "Creating new directory $EnrDirName"
    mkdir $EnrDirName
    cd $EnrDirName
    
    #######################################################
    # compute the reference solution
    
    # create directory for reference solution and go to it
    if [ -d  $RefSolnDirName ]; then
        echo "Deleting old directory $RefSolnDirName"
        rm -rf $RefSolnDirName
    fi
    echo "Creating new directory $RefSolnDirName"
    mkdir $RefSolnDirName
    cd $RefSolnDirName
    
    # copy template file
    cp $BaseFilePathReference ./
    echo "Copied template input file for reference solution into new directory."
    
    # modify the input file for the current run
    echo -n "Changing parameters ... "
    ReplaceString="s/i-0000/${EnrichmentValues[$e]}/1"
    sed -i $ReplaceString $ReferenceFileName.cpp
    ReplaceString="s/i-0002/${ReferenceRefinement}/1"
    sed -i $ReplaceString $ReferenceFileName.cpp
    echo "Done." 
    
    if [ "$UseDebug" -eq 1 ]; then
        # compile the input file
        echo -n "Compiling input file ... "
        $MORISROOT/share/scripts/create_shared_object.sh . build_dbg $ReferenceFileName >& log_compile.txt
        echo "Done."
        
        # run it
        echo -n "Perform simulation ... "
        $MRD ./$ReferenceFileName.so >& log.txt
        echo "Done."
    else
        # compile the input file
        echo -n "Compiling input file ... "
        $MORISROOT/share/scripts/create_shared_object.sh . build_opt $ReferenceFileName >& log_compile.txt
        echo "Done."
        
        # run it
        echo -n "Perform simulation ... "
        $MRO ./$ReferenceFileName.so >& log.txt
        echo "Done."
    fi
    
    # extract output values
    L2_REF="$(grep "L2_REF" log.txt | awk '{print $10}')"
    H1s_REF="$(grep "H1s_REF" log.txt | awk '{print $10}')"
    Num_DOFs="$(grep "Model - Total number of DOFs" log.txt | awk '{print $9}')"
    echo "$Num_DOFs,$L2_REF,$H1s_REF,$ReferenceElementSize" > ../errors.csv
    
    # go one back down in directory structure
    cd ..
    
    
    #######################################################
    # prepare and run coarser solutions 
    for i in ${!RefinementValues[@]}; do
        
        # create directory for refinement value
        RefineDirName="R_${RefinementValues[$i]}"
        echo ""
        echo "------------------------------------------------------"
        echo "Prepare runs for $RefineDirName"
        if [ -d  $RefineDirName ]; then
            echo "Deleting old directory $RefineDirName"
            rm -rf $RefineDirName 
        fi
        echo "Creating new directory $RefineDirName"
        mkdir $RefineDirName
        cd $RefineDirName
        
        #######################################################
        # compute the coarse solution
        
        # create directory for reference solution and go to it
        if [ -d  $CoarseSolnDirName ]; then
            echo "Deleting old directory $CoarseSolnDirName"
            rm -rf $CoarseSolnDirName
        fi
        echo "Creating new directory $CoarseSolnDirName"
        mkdir $CoarseSolnDirName
        cd $CoarseSolnDirName
        
        # copy template file
        cp $BaseFilePathCoarse ./
        echo "Copied template of input file for COARSE solution into new directory."
        
        # modify the input file for the current run
        echo -n "Changing parameters ... "
        ReplaceString="s/i-0000/${EnrichmentValues[$e]}/1"
        sed -i $ReplaceString $CoarseFileName.cpp
        ReplaceString="s/i-0001/${RefinementValues[$i]}/1"
        sed -i $ReplaceString $CoarseFileName.cpp
        ReplaceString="s/i-0002/${ReferenceRefinement}/1"
        sed -i $ReplaceString $CoarseFileName.cpp
        echo "Done." 
        
        if [ "$UseDebug" -eq 1 ]; then
            # compile the input file
            echo -n "Compiling input file ... "
            $MORISROOT/share/scripts/create_shared_object.sh . build_dbg $CoarseFileName >& log_compile.txt
            echo "Done."
            
            # run it
            echo -n "Perform simulation ... "
            $MRD ./$CoarseFileName.so >& log.txt
            echo "Done."
        else
            # compile the input file
            echo -n "Compiling input file ... "
            $MORISROOT/share/scripts/create_shared_object.sh . build_opt $CoarseFileName >& log_compile.txt
            echo "Done."
            
            # run it
            echo -n "Perform simulation ... "
            $MRO ./$CoarseFileName.so >& log.txt
            echo "Done."
        fi
        
        
        # go one back down in directory structure
        cd ..
        
        #######################################################
        # compare solutions
        
        # create directory for reference solution and go to it
        if [ -d  $CompareSolnDirName ]; then
            echo "Deleting old directory $CompareSolnDirName"
            rm -rf $CompareSolnDirName
        fi
        echo "Creating new directory $CompareSolnDirName"
        mkdir $CompareSolnDirName
        cd $CompareSolnDirName
        
        # copy template file
        cp $BaseFilePathCompare ./
        echo "Copied template of input file for COMPARE solution into new directory."
        
        # modify the input file for the current run
        echo -n "Changing parameters ... "
        ReplaceString="s/i-0000/${EnrichmentValues[$e]}/1"
        sed -i $ReplaceString $CompareFileName.cpp
        ReplaceString="s/i-0001/${RefinementValues[$i]}/1"
        sed -i $ReplaceString $CompareFileName.cpp
        ReplaceString="s/i-0002/${ReferenceRefinement}/1"
        sed -i $ReplaceString $CompareFileName.cpp
        echo "Done." 
        
        if [ "$UseDebug" -eq 1 ]; then
            # compile the input file
            echo -n "Compiling input file ... "
            $MORISROOT/share/scripts/create_shared_object.sh . build_dbg $CompareFileName >& log_compile.txt
            echo "Done."
            
            # run it
            echo -n "Perform simulation ... "
            $MRD ./$CompareFileName.so >& log.txt
            echo "Done."
        else
            # compile the input file
            echo -n "Compiling input file ... "
            $MORISROOT/share/scripts/create_shared_object.sh . build_opt $CompareFileName >& log_compile.txt
            echo "Done."
            
            # run it
            echo -n "Perform simulation ... "
            $MRO ./$CompareFileName.so >& log.txt
            echo "Done."
        fi

        
        # extract output values
        L2_Global="$(grep "L2_Global" log.txt | awk '{print $9}')"
        H1s_Global="$(grep "H1s_Global" log.txt | awk '{print $9}')"
        Num_DOFs="$(grep "Model - Total number of DOFs" log.txt | awk '{print $9}')"
        echo "$Num_DOFs,$L2_Global,$H1s_Global,${ElementSizes[$i]}" >> ../../errors.csv
        
        # go one back down in directory structure
        cd ..
        
        #######################################################
        
        # go one back down in directory structure
        cd ..
        
    # end: loop over ghost values
    done
    
    # go one back down in directory structure
    cd ..
    
# end: loop over enrichment values
done

#############################################################
echo ""
echo "======================================================"
echo "============= Done, all runs completed. =============="
echo "======================================================"
#############################################################
