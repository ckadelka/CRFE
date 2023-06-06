#!/bin/bash

#PBS -lwalltime=24:00:00
#PBS -lnodes=1:ppn=8
#PBS -W group_list=ithaca
#PBS -q normal_q

#add your modules
module purge
module load gcc python atlas

# This scripts loads multiple single processor tasks on one BlueRidge node.
# Each task uses input files pulled from a tarball. Other approaches would 
# be to specify different command line inputs for each task, or to simply 
# run the same task multiple times (e.g. for Monte Carlo simulations)

# To set up:
# Create tarballs with all the files associated with the task
# Name the tarballs with a counter(number) embedded in the name, 
# example, data_files_1.tar data_files_2.tar and so on

# $PBS_NP is the total number of processors allocated (16)
# $np is the counter for creating different directories for 
# different task with a name having $np in it

# To customize and run multiple tasks requiring more than one processor
# modify $PBS_NP appropriately

np=1
while [ $np -le $PBS_NP ]
do
  # Run from this directory
#rundir=$WORK/$np
  rundir=$PBS_O_WORKDIR

  # Change to the directory where the job was submitted
  cd $PBS_O_WORKDIR
  
  # Output the task number (1, 2, 3, etc)
  echo "Starting task #$np running from directory $rundir"

  # Create a new directory for the task and copy data files into it  
#  mkdir -p $rundir
#  cp data_files_$np.tar $rundir/

  # Enter the new task directory and untar the data files
#  cd $rundir
#  tar xf data_files_$np.tar
#  rm data_files_$np.tar

  # Launching the application with '&' at the end sends the process in background
  # Users should make sure to customize the name of the executable (here a.out)
  # Output will be written to 1.out for task 1, 2.out for task 2, etc
# ./a.out > $np.out &

nsim=20
w=10

if [ $np -eq 1 ]
then
python multiple_runs_of_MCMC_collect_data04.py -n$nsim -w$w -s0 -b10 &
elif [ $np -eq 2 ]
then
python multiple_runs_of_MCMC_collect_data04.py -n$nsim -w$w -s0 -b5 &
elif [ $np -eq 3 ]
then
python multiple_runs_of_MCMC_collect_data04.py -n$nsim -w$w -s0 -b2 &
else
python multiple_runs_of_MCMC_collect_data04.py -n$nsim -w$w -s1 -b1 &

fi



  # Increment task counter
  np=$((np+1))
done

# It is important to have wait at the end for all the processes to finish
wait

exit;
