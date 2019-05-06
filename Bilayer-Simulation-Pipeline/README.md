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

Here is a simple guideline on how to use our scripts:

1. Install the simulationo pipeline through calling the following command:
	> $ python3 setup.py
	* You can set the directory for installation through changing the parentDir entry in configuration.ini

2. Source GROMACS if necessary: usually this line of command works:
	> $ source /usr/local/gromacs/bin/GMXRC

3. To set up your simulation configuration, you do the following:
    * Directly edit the configuration.ini, 
    * Please strictly follow the format for the configuration file: you can add configuration for more simulations, but each simulation configuration should be uniform

4. To do simulations, please copy the following files into the folder where you want the simulation results to be held
	* configuration.ini
	* simulate.py or simulate-multi.py (the latter is for multi-processing multiple simulations at the same time

5. Directly run the following command for the simulation to run:
	*python3 simulate.py (or simulate-multi.py
