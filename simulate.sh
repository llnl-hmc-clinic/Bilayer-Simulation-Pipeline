#!/bin/bash

# command line arguments for insane.py 
args=$1

# info for folder naming and mdp file editing
descr=$2    # tc-grps
n=$3
asym=$4

# replacement strings for mdp files
s1='1.0 '   # tau_t
s2='320 '   # ref_t
s3='3e-4 '  # compressibility
s4='1.0 '   # ref_p

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

# energy minimization

mkdir "$descr $asym"
cd "$descr $asym"
mkdir em
cd em
../../files/insane.py $args
cat ../../files/header.txt top.top > topol.top
sed -i '' '3d' topol.top
cp ../../files/minimization.mdp ../em/
sed -i '' "s/REPLACE0/$descr/" minimization.mdp
sed -i '' "s/REPLACE1/$str1/" minimization.mdp
sed -i '' "s/REPLACE2/$str2/" minimization.mdp
sed -i '' "s/REPLACE3/$str3/" minimization.mdp
sed -i '' "s/REPLACE4/$str4/" minimization.mdp
gmx grompp -f minimization.mdp -c bilayer.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -v
cd ..

# 1fs

mkdir 1fs
cp em/em.gro em/topol.top 1fs
cd 1fs
cp ../../files/martini_v2.x_new-rf.1fs.mdp ../1fs/
sed -i '' "s/REPLACE0/$descr/" martini_v2.x_new-rf.1fs.mdp
sed -i '' "s/REPLACE1/$str1/" martini_v2.x_new-rf.1fs.mdp
sed -i '' "s/REPLACE2/$str2/" martini_v2.x_new-rf.1fs.mdp
sed -i '' "s/REPLACE3/$str3/" martini_v2.x_new-rf.1fs.mdp
sed -i '' "s/REPLACE4/$str4/" martini_v2.x_new-rf.1fs.mdp
gmx grompp -f martini_v2.x_new-rf.1fs.mdp -c em.gro -p topol.top -o 1fs.tpr
gmx mdrun -deffnm 1fs -v -nt 4 -dlb yes
cd ..

# 5fs 

mkdir 5fs
cp 1fs/1fs.gro 1fs/topol.top 5fs
cd 5fs
cp ../../files/martini_v2.x_new-rf.5fs.mdp ../5fs/
sed -i '' "s/REPLACE0/$descr/" martini_v2.x_new-rf.5fs.mdp
sed -i '' "s/REPLACE1/$str1/" martini_v2.x_new-rf.5fs.mdp
sed -i '' "s/REPLACE2/$str2/" martini_v2.x_new-rf.5fs.mdp
sed -i '' "s/REPLACE3/$str3/" martini_v2.x_new-rf.5fs.mdp
sed -i '' "s/REPLACE4/$str4/" martini_v2.x_new-rf.5fs.mdp
gmx grompp -f martini_v2.x_new-rf.5fs.mdp -c 1fs.gro -p topol.top -o 5fs.tpr
gmx mdrun -deffnm 5fs -v -nt 4 -dlb yes
cd ..

# 15fs

mkdir 15fs
cp 5fs/5fs.gro 5fs/topol.top 15fs
cd 15fs
cp ../../files/martini_v2.x_new-rf.15fs.mdp ../15fs/
sed -i '' "s/REPLACE0/$descr/" martini_v2.x_new-rf.15fs.mdp
sed -i '' "s/REPLACE1/$str1/" martini_v2.x_new-rf.15fs.mdp
sed -i '' "s/REPLACE2/$str2/" martini_v2.x_new-rf.15fs.mdp
sed -i '' "s/REPLACE3/$str3/" martini_v2.x_new-rf.15fs.mdp
sed -i '' "s/REPLACE4/$str4/" martini_v2.x_new-rf.15fs.mdp
gmx grompp -f martini_v2.x_new-rf.15fs.mdp -c 5fs.gro -p topol.top -o 15fs.tpr
gmx mdrun -deffnm 15fs -v -nt 4 -dlb yes
cd ..

# 20fs

mkdir 20fs
cp 15fs/15fs.gro 15fs/topol.top 20fs
cd 20fs
cp ../../files/martini_v2.x_new-rf.prod_run.mdp ../20fs/
sed -i '' "s/REPLACE0/$descr/" martini_v2.x_new-rf.prod_run.mdp
sed -i '' "s/REPLACE1/$str1/" martini_v2.x_new-rf.prod_run.mdp
sed -i '' "s/REPLACE2/$str2/" martini_v2.x_new-rf.prod_run.mdp
sed -i '' "s/REPLACE3/$str3/" martini_v2.x_new-rf.prod_run.mdp
sed -i '' "s/REPLACE4/$str4/" martini_v2.x_new-rf.prod_run.mdp
gmx grompp -f martini_v2.x_new-rf.prod_run.mdp -c 15fs.gro -p topol.top -o 20fs.tpr
gmx mdrun -deffnm 20fs -v -nt 4 -dlb yes

# indexing

mkdir analysis
cd analysis
(echo del 0-100; echo a PO4 ROH; echo name 0 headgroups; echo q) | gmx make_ndx -f ../20fs.gro -o index.ndx
gmx trjconv -f ../20fs.xtc -s ../20fs.tpr -o traj-center-mol.xtc -center -n index.ndx -pbc mol

# analysis

fatslim membranes -c ../20fs.gro -n index.ndx --output-index bilayer_leaflet.ndx
fatslim thickness -c ../20fs.gro -n index.ndx -t ../20fs.xtc --plot-thickness thickness.xvg
fatslim apl -c ../20fs.gro -n index.ndx -t ../20fs.xtc --plot-apl apl.xvg --plot-area area.xvg

# plots

mkdir plots
cd plots
../../../../files/plot.py ../thickness.xvg -o thickness.pdf
../../../../files/plot.py ../area.xvg -o area.pdf
../../../../files/plot.py ../apl.xvg -o apl.pdf
