"""

HMC Clinic Project with Lawrence Livermore National Lab 
Final Version 4/19/2019

This script runs a batch of bilayer simulations with a common compositional
leaflet asymmetry across a range of asymmetries. If the two leaflets do not
have the same composition, the script will also run symmetrical simulations
of the bilayer.

"""

import sys, subprocess

### environmental parameters ###

x, y, z = (12, 12, 10)  # box dimensions        
salt = 0.15             # molar salt concentration

### bilayer parameters ###

numLipids = 222 # number of lipids in Upper Leaflet
lipids = ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"]
upper = [243, 121, 20,  61, 242,   0,  0, 313]
lower = [139,  75, 54, 161, 108, 161, 22, 280]
#asymList = range(-30, -25)
asymList = [10]

### helper function that initializes arguments to pass to insane ###
def args_for_insane(upperCounts, lowerCounts):
	args = f" -pbc square -sol W -salt {salt} -x {x} -y {y} -z {z} -num {numLipids} -o bilayer.gro -p top.top "
	i = 0
	for l in lipids:
		if lower[i] > 0:
			args += f" -l {l}:{lowerCounts[i]} "
		if upper[i] > 0:
			args += f" -u {l}:{upperCounts[i]} "
		i += 1
	args += " -asym "
	return args

queue = [] # job queue

### add symmetrical runs ###

run_num = 0
if upper != lower: # only add symmetrical runs if upper/lower leaflets have different composition
	compositions = [(upper, lower), (upper, upper), (lower, lower)]
	for comp in compositions:
		args = args_for_insane(comp)
		for j in range(0, 4): # run number
			runDirectory = "symmetric/run" + str(run_num)
			argsForInsane = args + str(0)
			queue.append(argsForInsane)
			queue.append(runDirectory)
			run_num += 1

### add asymmetrical runs ###

run_num = 0
args = args_for_insane(upper, lower)
for asym in asymList:
	for j in range(0, 4):
		runDirectory = "run" + str(run_num)
		argsForInsane = args + str(asym)
		queue.append(argsForInsane)
		queue.append(runDirectory)
		run_num += 1 

### run simulations ###

while len(queue) > 0:
	command = ["sbatch", "setup.sh"] + queue[:8]
	subprocess.call(command)
	queue = queue[8:]

sys.exit(0)
