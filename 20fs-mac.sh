#!/bin/sh

descr=$1
counts=$2
runnum=$3
n=$4
str5=$5

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
	sed -i " " "s/REPLACE0/$descr/" $1
	sed -i " " "s/REPLACE1/$str1/" $1
	sed -i " " "s/REPLACE2/$str2/" $1
	sed -i " " "s/REPLACE3/$str3/" $1
	sed -i " " "s/REPLACE4/$str4/" $1
	sed -i " " "s/REPLACE5/$str5/" $1

}


### 20fs ###
cd $runnum
mkdir 20fs
cp 15fs/15fs.gro 15fs/topol.top 20fs
cd 20fs
cp ../../files/martini_v2.x_new-rf.prod_run.mdp ../20fs/
make_mdp "martini_v2.x_new-rf.prod_run.mdp"
gmx grompp -f martini_v2.x_new-rf.prod_run.mdp -c 15fs.gro -p topol.top -o 20fs.tpr -maxwarn 1
gmx mdrun -deffnm 20fs -v -nt 8
