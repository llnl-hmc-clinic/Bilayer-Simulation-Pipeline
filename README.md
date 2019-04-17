	Welcome to the usage of our lipid bilayer simulation script! 

This script serves the purpose of doing multiple simulations in one run with simplified commands. Also, you can do analysis of lipid bilayer properties with our script.

To use the script, there are some packages and softwares that are required:

1. python3.6 both installed

2. GROMACS
    http://www.gromacs.org/

3. INSert membrANE (insane.py)
It has been included in the package, so usually you don't need to worry about it.

For analysis purposes, you should install also fatslim
http://fatslim.github.io/


To use our script, you need to:

1. Activate GROMACS if necessary: usually this line of command works:
	$source /usr/local/gromacs/bin/GMXRC

2. To set up your simulation configuration, you can choose one of the following methods:
    * Directly edit the configuration.csv, 
    * Make changes to mkcsv.py and then run this script to generate your csv file.
    * You can also use mkcsv.py to generate a template for your simulation configuration in case your configuration file is deleted by accident.
    * Please strictly follow the format for the configuration file: you can add configuration for more simulations, but each simulation configuration should be uniform

3. To do simulations:
	If you have no need for multiprocessing, run 
		> $ python3 simulate.py

	If you are multiprocessing the simulations, run
		> $ python3 simulate-multi.py

	If you are multiprocessing on sbatch, run
		> $ python3 simulate-super.py

Our pipeline is for simulations from energy minimization to your production run.