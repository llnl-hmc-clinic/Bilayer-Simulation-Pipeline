#!/bin/bash

#################################################################################
###### HMC Clinic Project, Final Version 4/19/2019                         ######
######                                                                     ######
###### This sbatch script runs the long time-scale simulation used for     ######
###### analysis. It takes as its first argument the number of recursive    ######
###### calls to make, given the 16 hour limit for jobs on syrah, to ensure ######
###### the full time course is simulated. The next four arguments specify  ######
###### the working directories where the setup script has been run. Once   ######
###### the simulation is complete, we call the analysis script and exit.   ######
######                                                                     ######
#################################################################################

### Use: sbatch run-syrah.sh X R1 R2
### where X is the number of times you want it to repeat.

#SBATCH -J HM_init
#SBATCH -A bbs
#SBATCH -N 1
#SBATCH -t 16:00:00
#SBATCH -p pbatch
#SBATCH --export=ALL

echo "Job start running"

# pbatch - pdebug
# 16:00:00 - 0.1 - 16

# Set flags for mdrun (-cpi is set by the script when needed, no need to to include it here)
FLAGS='-dlb yes -maxh 16 -deffnm 20fs -cpi '

export g_gmx="/p/lustre1/muntere/gromacs-5.1.4/bin/gmx"

# Decrease maximum number of iterations ($maxrun arg passed to script)
if [ -z "$1" ] 
then  
  maxrun=0
else
  maxrun=$1
fi
let maxrun--

echo "Run jobs $2 $3 $4 $5 for $1 times"

# Set working directory
echo "The jobs are ran here:"
pwd

export GROMACS_HOME=/p/lustre1/muntere/gromacs-5.1.4
export gromacs="${GROMACS_HOME}/bin/gmx_mpi mdrun"

source $GROMACS_HOME/bin/GMXRC.bash

echo "Running job $SLURM_JOB_ID with #nodes 1 and #prog 4 maxrun $maxrun"

# Start next job, but on hold waiting for this one
if [ $maxrun -gt 0 ] 
then
    echo "sbatch --dependency=afterok:$SLURM_JOB_ID run-syrah.sh $maxrun $2 $3 $4 $5 >> job.id"
    sbatch --dependency=afterok:$SLURM_JOB_ID run-syrah.sh $maxrun $2 $3 $4 $5 >> job.id
fi

cd "$2/20fs"
  echo "Running jobs $2"
  $g_gmx mdrun -deffnm 20fs -nt 4 -dlb yes -maxh 16 -cpi >> mdrun.log 2>&1 &
cd ../../
cd "$3/20fs"
  echo "Running jobs $3"
  $g_gmx mdrun -deffnm 20fs -nt 4 -dlb yes -maxh 16 -cpi >> mdrun.log 2>&1 &
cd ../../
cd "$4/20fs"
  echo "Running jobs $4"
  $g_gmx mdrun -deffnm 20fs -nt 4 -dlb yes -maxh 16 -cpi >> mdrun.log 2>&1 &
cd ../../
cd "$5/20fs"
  echo "Running jobs $5"
  $g_gmx mdrun -deffnm 20fs -nt 4 -dlb yes -maxh 16 -cpi >> mdrun.log 2>&1 &
cd ../../

wait
sbatch analysis.sh $2 $3 $4 $5


