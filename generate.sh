#!/bin/bash

### variables ###

args=$1
descr=$2
n=$3
asym=$4
counts=$5
runnum=$6
Utextnote=$7
Ltextnote=$8
Atextnote=$9

### energy minimization ###

mkdir $runnum
cd $runnum
echo $Utextnote >description.txt
echo $Ltextnote >> description.txt
echo $Atextnote >> description.txt
mkdir em
cd em
python2.7 ../../files/insane.py $args
cat ../../files/header.txt top.top > topol.top
sed -i '5d' topol.top
cp ../../files/minimization.mdp ../em/
(echo del 1-200; echo "r W | r NA+ | r CL-"; echo name 1 Solvent; echo !1; echo name 2 Membrane; echo q) | gmx make_ndx -f bilayer.gro -o index.ndx
gmx grompp -f minimization.mdp -c bilayer.gro -p topol.top -n index.ndx -o em.tpr
gmx mdrun -deffnm em -v -nt 2 -dlb yes
cd ..

### 1fs ###

mkdir 1fs
cp em/em.gro em/topol.top 1fs
cd 1fs
cp ../../files/martini_v2.x_new-rf.1fs.mdp ../em/index.ndx ../1fs/
gmx grompp -f martini_v2.x_new-rf.1fs.mdp -c em.gro -p topol.top -n index.ndx -o 1fs.tpr
gmx mdrun -deffnm 1fs -v -nt 4 -dlb yes
cd ..

### 5fs ### 

mkdir 5fs
cp 1fs/1fs.gro 1fs/topol.top 5fs
cd 5fs
cp ../../files/martini_v2.x_new-rf.5fs.mdp ../1fs/index.ndx ../5fs/
gmx grompp -f martini_v2.x_new-rf.5fs.mdp -c 1fs.gro -p topol.top -n index.ndx -o 5fs.tpr
gmx mdrun -deffnm 5fs -v -nt 6 -dlb yes
cd ..

### 15fs ###

mkdir 15fs
cp 5fs/5fs.gro 5fs/topol.top 15fs
cd 15fs
cp ../../files/martini_v2.x_new-rf.15fs.mdp ../5fs/index.ndx ../15fs/
gmx grompp -f martini_v2.x_new-rf.15fs.mdp -c 5fs.gro -p topol.top -n index.ndx -o 15fs.tpr
gmx mdrun -deffnm 15fs -v -nt 6 -dlb yes
cd ..

