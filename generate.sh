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

### helper function for .mdp file generation ###

s1='1.0 '
s2='320 '
s3='3e-4 '
s4='1.0 '

str1=''
str2=''
str3=''
str4=''
while [ $n -ge 1 ]
do
	str1+="$s1"
	str2+="$s2"
	str3+="$s3"
	str4+="$s4"
	((n--))
done


function make_mdp {
	sed -i "s/REPLACE0/$descr/" $1
	sed -i "s/REPLACE1/$str1/" $1
	sed -i "s/REPLACE2/$str2/" $1
	sed -i "s/REPLACE3/$str3/" $1
	sed -i "s/REPLACE4/$str4/" $1
}

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
sed -i 5d topol.top
cp ../../files/minimization.mdp ../em/
make_mdp "minimization.mdp"
gmx grompp -f minimization.mdp -c bilayer.gro -p topol.top -o em.tpr -maxwarn 1
gmx mdrun -deffnm em -v
cd ..

### 1fs ###

mkdir 1fs
cp em/em.gro em/topol.top 1fs
cd 1fs
cp ../../files/martini_v2.x_new-rf.1fs.mdp ../1fs/
make_mdp "martini_v2.x_new-rf.1fs.mdp"
gmx grompp -f martini_v2.x_new-rf.1fs.mdp -c em.gro -p topol.top -o 1fs.tpr -maxwarn 1
gmx mdrun -deffnm 1fs -v -nt 2 -dlb yes
cd ..

### 5fs ### 

mkdir 5fs
cp 1fs/1fs.gro 1fs/topol.top 5fs
cd 5fs
cp ../../files/martini_v2.x_new-rf.5fs.mdp ../5fs/
make_mdp "martini_v2.x_new-rf.5fs.mdp"
gmx grompp -f martini_v2.x_new-rf.5fs.mdp -c 1fs.gro -p topol.top -o 5fs.tpr -maxwarn 1
gmx mdrun -deffnm 5fs -v -nt 2 -dlb yes
cd ..

### 15fs ###

mkdir 15fs
cp 5fs/5fs.gro 5fs/topol.top 15fs
cd 15fs
cp ../../files/martini_v2.x_new-rf.15fs.mdp ../15fs/
make_mdp "martini_v2.x_new-rf.15fs.mdp"
gmx grompp -f martini_v2.x_new-rf.15fs.mdp -c 5fs.gro -p topol.top -o 15fs.tpr -maxwarn 1
gmx mdrun -deffnm 15fs -v -nt 2 -dlb yes
cd ../..

