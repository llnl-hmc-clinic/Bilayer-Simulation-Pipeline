	Welcome to the usage of our lipid bilayer simulation script! 

This script serves the purpose of doing multiple simulations in one run with simplified commands. Also, you can do analysis of lipid bilayer properties with our script.

To use the script, there are some packages and softwares that are required:

1: python 2.7 and python3.6 both installed

2: GROMACS
http://www.gromacs.org/

3: INSert membrANE (insane.py)
It has been included in the package, so usually you don't need to worry about it.

For analysis purposes, you can choose to use fatslim, which is our default analysis tool for some property analysis

http://fatslim.github.io/


To use our script, you need to:

0: (If you use this script for the first time) Call the following command on terminal: 
	$chmod u+x simulate.py generate.sh generate-mac.sh 20fs.sh 20fs-mac.sh

1: Activate GROMACS if necessary: usually this line of command works:
	$source /usr/local/gromacs/bin/GMXRC

2: To set up your simulation configuration, you can either:
	Directly edit the configuration.csv, or
	Make changes to mkcsv.py and then run this script to generate your csv file.
	You can also use mkcsv.py to generate a template for your code in case your configuration file is deleted by accident.
	Please strictly follow the format for the configuration file: each simulation needs 5 entries, and there should not be a blank line between different simulations.

3: To do simulations:
	If you want to do relaxation of the membrane, which is doing simulation of energy minimization, 1fs, 5fs and 15fs timesteps, run:
	$python3 simulate.py relax
	If you want to only run the 20fs simulation on the basis of previous membrane relaxation, run:
	$python3 simulate.py 20fs
	If you want to do a complete simulation, run:
	$python3 simulate.py
(this is for simulations from energy minimization to 15fs runs)

4: To analyze the simulation results: 
	$python3 analyze.py
