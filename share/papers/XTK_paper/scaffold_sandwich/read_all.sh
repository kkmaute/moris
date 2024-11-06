#!/bin/bash

### if shell script not accepted by unix, do:
#   chmod +x read_all.sh
#   sed -i -e 's/\r$//' read_all.sh

### to run script:
#   ./read_all.sh

#############################################################
echo ""
echo "================================================"
echo "======== Read results from all files... ========"
echo "================================================"
echo ""
#############################################################

BaseDirectory=$PWD
echo "Base directory = $BaseDirectory "

Studies=("run_1" "run_2" "run_3")
SubStudies=("size" "weak" "strong_large" "strong_small")

for i in ${!Studies[@]}; do

    cd ${Studies[$i]}

    for j in ${!SubStudies[@]}; do
    
        DirName="${SubStudies[$j]}_scaling"
        ScriptName="read_results_${SubStudies[$j]}.sh"
        cd $DirName
        chmod +x $ScriptName
        sed -i -e 's/\r$//' $ScriptName
        ./$ScriptName
        cd ..

    done
    
    cd ..

done