#!/bin/bash

### if shell script not accepted by unix, do:
#   chmod +x read_results_strong.sh
#   sed -i -e 's/\r$//' read_results_strong.sh

### to run script:
#   ./read_results_strong.sh

#############################################################
echo ""
echo "================================================"
echo "======== Read results from log files... ========"
echo "================================================"
echo ""
#############################################################

BaseDirectory=$PWD
echo "Base directory: $BaseDirectory "
export LogFileName="log.txt"
echo "File name for log files: $LogFileName "

# create file for timing results
export TimingFileName="timings.csv"
echo "DecompTime,EnrichTime,GhostTime" > "$TimingFileName"

export MemoryFileName="memory.csv"
echo "DecompMemIn,EnrichMemIn,GhostMemIn,GhostMemOut" > "$MemoryFileName"

### For Alpine:
NumProcs=("1" "2" "4" "8" "16" "32" "64")


for i in ${!NumProcs[@]}; do
    
    # make sure to start from the base directory
    cd $BaseDirectory

    # go to directory
    NewDirName="np_${NumProcs[$i]}"
    cd $NewDirName

    # output
    echo ""
    echo "======================================================"
    echo "Read data from $NewDirName ... "
    
    # create output array
    TimingVals=()
    MemoryVals=()
    
    # find decomp in memory consumption
    line_number=$(grep -n "|__XTK - Overall - Run " "$LogFileName" | cut -d: -f1)
    target_line=$((line_number + 1))
    line=$(sed -n "${target_line}p" "$LogFileName")
    DecompMemIn=$(echo "$line" | awk '
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[+-]?[0-9]*[.]?[0-9]+$/) {
                    print $i
                    exit
                }
            }
        }
    ')
    
    # find the decomp time
    line_number=$(grep -n "|__XTK - NoType - Enrichment" "$LogFileName" | cut -d: -f1)
    target_line=$((line_number - 9))
    line=$(sed -n "${target_line}p" "$LogFileName")
    DecompTime=$(echo "$line" | awk '
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[+-]?[0-9]*[.]?[0-9]+$/) {
                    print $i
                    exit
                }
            }
        }
    ')
    target_line=$((line_number - 3))
    line=$(sed -n "${target_line}p" "$LogFileName")
    EnrichMemIn=$(echo "$line" | awk '
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[+-]?[0-9]*[.]?[0-9]+$/) {
                    print $i
                    exit
                }
            }
        }
    ')
    
    # find the enrichment time / enrichment out mem
    line_number=$(grep -n "|__XTK - NoType - Create high to low double side sets" "$LogFileName" | cut -d: -f1)
    target_line=$((line_number - 5))
    line=$(sed -n "${target_line}p" "$LogFileName")
    EnrichTime=$(echo "$line" | awk '
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[+-]?[0-9]*[.]?[0-9]+$/) {
                    print $i
                    exit
                }
            }
        }
    ')
    target_line=$((line_number - 3))
    line=$(sed -n "${target_line}p" "$LogFileName")
    GhostMemIn=$(echo "$line" | awk '
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[+-]?[0-9]*[.]?[0-9]+$/) {
                    print $i
                    exit
                }
            }
        }
    ')
    
    # find the ghost time / ghost out memory consumption
    line_number=$(grep -n "|_ - Registered mesh pair #0" "$LogFileName" | cut -d: -f1)
    target_line=$((line_number - 4))
    line=$(sed -n "${target_line}p" "$LogFileName")
    GhostTime=$(echo "$line" | awk '
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[+-]?[0-9]*[.]?[0-9]+$/) {
                    print $i
                    exit
                }
            }
        }
    ')
    target_line=$((line_number - 2))
    line=$(sed -n "${target_line}p" "$LogFileName")
    GhostMemOut=$(echo "$line" | awk '
        {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^[+-]?[0-9]*[.]?[0-9]+$/) {
                    print $i
                    exit
                }
            }
        }
    ')
    
    
    # write timing results to .csv
    cd ..
    echo "$DecompTime,$EnrichTime,$GhostTime" >> "$TimingFileName"
    echo "$DecompMemIn,$EnrichMemIn,$GhostMemIn,$GhostMemOut" >> "$MemoryFileName"

done


#############################################################
echo ""
echo "========================================="
echo "======== Done. All results read. ========"
echo "========================================="
#############################################################
