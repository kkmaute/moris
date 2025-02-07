#!/bin/bash

### if shell script not accepted by unix, do:
#   chmod +x run_brick_wall_study.sh
#   sed -i -e 's/\r$//' run_brick_wall_study.sh

### to run script:
#   ./run_brick_wall_study.sh >& log_run.txt &
#   tail -f log_run.txt

#############################################################
echo ""
echo "=============================================="
echo "======== Performing Brick Wall Study  ========"
echo "=============================================="
echo ""
#############################################################

echo "PID = $$"
BaseDirectory=$PWD
echo "Base directory = $BaseDirectory "
export FileName="brick_wall"
BaseFilePath="${BaseDirectory}/${FileName}.cpp"
export LogFileName="log.txt"

# Ghost Cases - identifier "i-00G"
GhostValues=("true" "false")

# Radius Cases - identifier "i-00R"
RadiusValues=("0.000" "0.001" "0.005" "0.010" "0.050" "0.100" "0.500" "1.000" "5.000")

# Off-Set Cases - identifiers "i-00X", "i-00Y", and "i-00Z"
OffSetValuesX=("0.0" "+6.1930e-03" "-1.0406e-01" "+6.4482e-02" "+2.2878e-01" "-1.4721e-02" "+3.7374e-03" "+2.1739e-01" "+6.5438e-03" "+3.8800e-01" "+6.0887e-01" "-5.9816e-02" "+9.0132e-01" "+6.8462e-01" "-1.4953e-02" "-9.4907e-01" "+6.3695e-01" "-1.3634e-01" "+3.2157e-01" "-3.5062e-02" "-2.0488e-01" "+7.4736e-01" "-1.8788e-01" "+3.1178e-01" "-3.6697e-01" "-8.9343e-02" "+7.6223e-02" "+1.9052e-01" "+2.0624e-01" "+1.5061e-01" "+6.1475e-02" "+6.4708e-02" "+9.1111e-01" "-2.8819e-01" "-7.3657e-01" "+3.2702e-02" "-3.6523e-02" "-2.1025e-01" "+5.5980e-02" "+9.6560e-02" "-2.5431e-01" "-5.6470e-01" "-6.5020e-01" "-4.6722e-01" "-1.1698e-03" "-3.0705e-01" "+4.2572e-01" "+1.1432e-03" "-1.7573e-02" "-1.7272e-01" "-3.4469e-01" "-7.4006e-02" "-1.5541e-01" "-4.9080e-03" "-7.3299e-01" "-2.3866e-01" "+2.2805e-01" "+8.5332e-01" "+5.2604e-02" "+1.1353e-02" "+7.1413e-01" "-3.2685e-02" "+5.0939e-01" "+2.0994e-02" "-7.6813e-01" "+3.4265e-01" "-6.9250e-02" "+1.5452e-02" "+3.2764e-01" "-3.1148e-01" "-9.7988e-03" "-1.4957e-01" )
OffSetValuesY=("0.0" "-6.1824e-01" "+9.4836e-01" "-9.2890e-02" "-9.3856e-01" "-1.8933e-01" "+9.0389e-02" "-2.0655e-01" "+2.7188e-02" "-8.8885e-01" "+8.4210e-02" "+1.8121e-03" "+8.3139e-01" "+1.6088e-02" "-5.1647e-02" "+6.4861e-03" "-3.5472e-02" "-8.3470e-01" "-4.0532e-01" "-4.1785e-02" "-2.6937e-02" "+1.6763e-01" "-2.1526e-02" "-2.4614e-02" "+4.4193e-02" "+9.9532e-01" "+3.5025e-01" "-8.2871e-02" "-5.8056e-01" "+4.2822e-01" "+6.2958e-01" "+7.1500e-03" "-1.7828e-01" "-6.6076e-01" "+1.0031e-01" "+7.2827e-02" "-4.2807e-03" "-3.3865e-01" "-1.5415e-01" "-9.5670e-01" "+2.4087e-01" "+2.9255e-01" "+4.2749e-02" "+2.9026e-01" "+4.1828e-01" "+7.3331e-02" "+3.0592e-01" "-2.5732e-02" "+7.3827e-01" "-1.0336e-02" "+4.6472e-02" "-2.9551e-03" "-2.1887e-01" "+7.7070e-01" "-8.9978e-01" "+5.5016e-03" "+8.5792e-01" "-1.8337e-01" "+1.9751e-01" "-3.5641e-02" "+2.2162e-01" "+4.9990e-01" "-3.3783e-01" "-1.7929e-01" "-3.7276e-02" "-5.5698e-03" "-3.0988e-01" "+1.1676e-01" "-1.3909e-01" "+1.2544e-02" "-8.5907e-02" "+2.6978e-03" )
OffSetValuesZ=("0.0" "-2.6696e-01" "-8.5719e-02" "+3.0471e-02" "+3.4686e-01" "-7.2004e-01" "+3.1572e-01" "+6.0832e-02" "-9.0486e-02" "-1.1127e-01" "+1.4921e-01" "+1.2508e-02" "+3.1607e-02" "-7.2657e-01" "-2.5228e-02" "-1.5440e-03" "-7.7290e-02" "+4.4108e-02" "-2.1930e-01" "+4.5659e-01" "+5.5548e-01" "-2.2734e-03" "+3.1529e-02" "+3.7184e-03" "-8.9810e-01" "+6.3343e-02" "-1.2337e-01" "+1.9259e-01" "+2.8323e-03" "+2.6664e-02" "-1.4457e-01" "-5.9767e-02" "+8.6350e-01" "-5.1949e-01" "-8.4085e-01" "+2.9831e-02" "+5.4279e-01" "+1.0073e-02" "-8.8807e-01" "-6.9202e-01" "+1.1173e-01" "+3.0794e-01" "-5.5991e-01" "-1.2423e-02" "-4.5745e-02" "+1.8188e-01" "+2.2276e-01" "-1.7319e-02" "+7.1798e-01" "+2.3162e-01" "-3.8878e-02" "-6.3408e-02" "+3.6910e-02" "-7.8625e-03" "+2.9118e-01" "+6.7039e-02" "-7.0875e-03" "-6.9787e-01" "-3.8879e-02" "+7.7991e-01" "-3.1146e-01" "-1.0444e-02" "+1.8132e-01" "+1.9258e-01" "-1.8108e-01" "+2.3257e-01" "-2.4368e-02" "-3.2982e-02" "-4.1293e-02" "-7.3928e-01" "-2.3495e-01" "+5.0868e-03" )


# ================================================================ #
for i in ${!GhostValues[@]}; do
    
    # make sure to start from the base directory
    cd $BaseDirectory

    # get the new folder name
    NewDirName="G_${GhostValues[$i]}"
    
    # output
    echo ""
    echo "===================================="
    echo "Prepare runs for Ghost: $NewDirName"
    
    # delete the directory if it already exists
    if [ -d  $NewDirName ]; then
        echo "Deleting old directory ${GhostValues[$i]}"
        rm -rf $NewDirName 
    fi
    
    # create directory for run
    echo "Creating new directory $NewDirName"
    mkdir $NewDirName
    
    # go to the new directory
    cd $NewDirName
    
    # ================================================================ #
    # process sub-directories 
    for j in ${!RadiusValues[@]}; do
    
        # get the new sub-folder name
        NewSubDirName="R_${RadiusValues[$j]}"
        
        # output
        echo ""
        echo "-------------------------------------------"
        echo "Prepare runs for Radius: ${RadiusValues[$j]}"
        
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
        
        # modify the input base parameters
        echo -n "Changing parameters ... "
        ReplaceString="s/i-00G/${GhostValues[$i]}/1"
        sed -i $ReplaceString $FileName.cpp
        ReplaceString="s/i-00R/${RadiusValues[$j]}/1"
        sed -i $ReplaceString $FileName.cpp
        echo "Done." 
        
        # initialize csv carrying values
        export ValueFileName="out.csv"
        echo "ghost,radius,off_x,off_y,off_z,lambda_max,lambda_min" > "$ValueFileName"
        count=0
        
        # ================================================================ #
        # run the various offset casess
        for k in ${!OffSetValuesX[@]}; do
        
            # create a temporary run  directory
            echo "---"
            echo "Iteration: $count"
            (( count++ ))
            echo -n "Prepare next offset value ... "
            RunDirName="run"
            if [ -d  $RunDirName ]; then
                rm -rf $RunDirName 
            fi
            mkdir $RunDirName
            cd $RunDirName
            cp ../$FileName.cpp ./
            
            # set the offset values
            ReplaceString="s/i-00X/${OffSetValuesX[$k]}/1"
            sed -i $ReplaceString $FileName.cpp
            ReplaceString="s/i-00Y/${OffSetValuesY[$k]}/1"
            sed -i $ReplaceString $FileName.cpp
            ReplaceString="s/i-00Z/${OffSetValuesZ[$k]}/1"
            sed -i $ReplaceString $FileName.cpp
            echo "Done."
        
            # compile the input file
            echo -n "Compiling input file ... "
            $MORISROOT/share/scripts/create_shared_object.sh . build_opt $FileName >& log_compile.txt
            echo "Done."

            # run it      
            echo -n "Perform simulation ... "
            mpirun -n 8 $MRO ./$FileName.so >& $LogFileName
            echo "Done."
            
            # extract values from log file file
            echo -n "Extract values from log file ... "
            OutStart="${GhostValues[$i]},${RadiusValues[$j]},${OffSetValuesX[$k]},${OffSetValuesY[$k]},${OffSetValuesZ[$k]}"

            # Max eigenvalue, + 4 lines after "|_SLEPc - KSP Solver: preonly", second numeric value (which has exponential format)
            line_number=$(grep -n "|_SLEPc - KSP Solver: preonly" "$LogFileName" | cut -d: -f1)
            target_line=$((line_number + 4))
            line=$(sed -n "${target_line}p" "$LogFileName")
            #echo $line
            LambdaMax=$( echo "$line" | awk '{print $13}')

            # Min eigenvalue, + 5 lines after "|_SLEPc - KSP Solver: fgmres", second numeric value (which has exponential format)
            line_number=$(grep -n "|_SLEPc - KSP Solver: fgmres" "$LogFileName" | cut -d: -f1)
            target_line=$((line_number + 5))
            line=$(sed -n "${target_line}p" "$LogFileName")
            #echo $line
            LambdaMin=$( echo "$line" | awk '{print $13}')

            echo "Done."
            
            # leave run directory and delete it
            cd ..
            rm -r $RunDirName
            
            # add values to .csv file
            echo "$OutStart,$LambdaMax,$LambdaMin" >> "$ValueFileName"
            
        # end: loop over offset values
        done
        # ================================================================ #
        
        # leave the radius directory
        cd ..
    
    # end: loop over radius values
    done
    # ================================================================ #
    
    # leave the ghost directory
    cd ..
    
# end: loop over ghost values
done
# ================================================================ #

#############################################################
echo ""
echo "======================================================"
echo "============= Done, all runs completed. =============="
echo "======================================================"
#############################################################