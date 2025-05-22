#!/bin/tcsh

### -----------------------------------------------------------------------
### P B S   D I R E C T I V E S  
### -----------------------------------------------------------------------

### This is the name which will be listed when looking at the list of jobs using squeue
### Note that only the first 8 characters of the name will be shown
#SBATCH --job-name=run_1
#SBATCH --output=size_scaling.%j.out

### Flags to specify that the whole node gets blocked for the job
#SBATCH --nodes=1 
#SBATCH --exclusive
#SBATCH --mem=750G
#SBATCH --cpu-freq=2000000
#SBATCH --ntasks-per-core=1
#SBATCH --nodelist=c3mem-a4-u34-[1-4],c3mem-a4-u36-[1-4],c3mem-a9-u7-[2,4]

### Time after which the script (and whatever was executed through it) is stopped
### if this flag is not specified your job may run forever if your running jobs on our own nodes
#SBATCH --time=2:30:00

### The below flags are used to define the allocation to run on
#SBATCH --partition=amem
#SBATCH --qos=mem

### slurm will send emails whenever a job starts, gets put on hold, finishes, or chrashes 
### to the email you specify here
#SBATCH --mail-user=niwu6047@colorado.edu
#SBATCH --mail-type=ALL

### -----------------------------------------------------------------------
### J O B   I N F O R M A T I O N  
### -----------------------------------------------------------------------

echo "PWD"=$PWD
echo "SLURM_JOBID"=$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo ""

### -----------------------------------------------------------------------
### SOURCE ENVIRONMENT
### -----------------------------------------------------------------------

### read all the environment variables defined for MORIS
source ~/.cshrc_moris_alpine

### -----------------------------------------------------------------------
### B E G I N N I N G   O F   E X E C U T I O N 
### -----------------------------------------------------------------------

###! IMPORTANT: run compile scripts for each scaling case, then create copies of the pre-compiled run-directory

### set name of the run-directory
RunDir="run_1"

### SIZE SCALING
echo -n "Run size scaling ... "
cd "$RunDir/size_scaling"
chmod +x run_size_scaling.sh
sed -i -e 's/\r$//' run_size_scaling.sh
./run_size_scaling.sh >& log_run_size_scaling.txt
cd ../..
echo "Done."

### STRONG SCALING LARGE
echo -n "Run strong scaling small ... "
cd "$RunDir/strong_scaling_large"
chmod +x run_strong_scaling_large.sh
sed -i -e 's/\r$//' run_strong_scaling_large.sh
./run_strong_scaling_large.sh >& log_run_strong_scaling_large.txt
cd ../..
echo "Done."

### STRONG SCALING SMALL
echo -n "Run strong scaling small ... "
cd "$RunDir/strong_scaling_small"
chmod +x run_strong_scaling_small.sh
sed -i -e 's/\r$//' run_strong_scaling_small.sh
./run_strong_scaling_small.sh >& log_run_strong_scaling_small.txt
cd ../..
echo "Done."

### WEAK SCALING
echo -n "Run weak scaling ... "
cd "$RunDir/weak_scaling"
chmod +x run_weak_scaling.sh
sed -i -e 's/\r$//' run_weak_scaling.sh
./run_weak_scaling.sh >& log_run_weak_scaling.txt
cd ../..
echo "Done."
