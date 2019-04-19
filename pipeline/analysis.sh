#!/bin/bash

#################################################################################
###### HMC Clinic Project, Final Version 4/19/2019                         ######
######                                                                     ######
###### This sbatch script performs analysis on up to four jobs at once.    ######
###### It takes four arguments for the working directories to run within,  ######
###### assuming they have properly completed set up and production run.    ######
######                                                                     ######
#################################################################################

### Use: sbatch analysis.sh run1 run2 run3 run4

#SBATCH -J HM_init
#SBATCH -A bbs
#SBATCH -N 1
#SBATCH -t 6:00:00
#SBATCH -p pbatch
#SBATCH --export=ALL

echo "Job start running"

# pbatch - pdebug
# 2:00:00 - 0.1 - 16

source '~/.bashrc'
source '/p/lustre1/muntere/gromacs-5.1.4/bin/gmx'

echo "Run job $1"

# Set working directory
echo "The jobs are ran here:"
pwd

export GROMACS_HOME=/p/lustre1/muntere/gromacs-5.1.4
export gromacs="${GROMACS_HOME}/bin/gmx_mpi mdrun"

source $GROMACS_HOME/bin/GMXRC.bash

flags () {
	return "-r $1/20fs/20fs.trr -f $1/20fs/20fs.xtc -s $1/20fs/20fs.tpr -n $1/20fs/index.ndx -g $1/20fs/20fs.gro -d $1/20fs/20fs.edr -m $1/20fs/mdout.mdp -p $1/20fs/topol.tpr -o $1/analysis/"
}

echo "Running jobs $1 $2 $3 $4"
python3 scripts/analysis.py $((flags $1)) &
python3 scripts/analysis.py $((flags $2)) &
python3 scripts/analysis.py $((flags $3)) &
python3 scripts/analysis.py $((flags $4)) &

wait


