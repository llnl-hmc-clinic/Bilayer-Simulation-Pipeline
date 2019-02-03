#!/bin/bash

### variables ###

descr=$1
counts=$2
runnum=$3
n=$4
str5=$5

### helper function for .mdp file generation ###

function make_mdp {
	sed -i "s/REPLACE5/$str5/" $1
}

### 20fs ###
cd $runnum
mkdir 20fs
cp 15fs/15fs.gro 15fs/topol.top 15fs/index.ndx 20fs
cd 20fs
cp ../../files/martini_v2.x_new-rf.prod_run.mdp ../20fs/
make_mdp "martini_v2.x_new-rf.prod_run.mdp"
gmx grompp -f martini_v2.x_new-rf.prod_run.mdp -c 15fs.gro -p topol.top -n index.ndx -o 20fs.tpr
gmx mdrun -deffnm 20fs -v -nt 8 -dlb yes
cd ../..
