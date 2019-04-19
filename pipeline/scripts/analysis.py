"""

HMC Clinic Project with Lawrence Livermore National Lab 
Final Version 4/19/2019

This script extends the mdreader class for lipid bilayers and includes analysis functions 
for computing thickness, area per lipid, area compressibility, order parameters, leaflet 
surface tension, diffusion coefficients, and cholesterol flip flop rate. Other functions 
included are meant to integrate the functionality of gromacs, fatslim and mdreader to allow
easy extension to other property analysis.

"""

import os, subprocess, sys, re, getopt
import math, statistics, shlex, mdreader
import numpy as np
from numpy import linspace, exp
from numpy.random import randn
from scipy import integrate
from scipy.interpolate import UnivariateSpline
from scipy.optimize import leastsq

usage = """
usage: $ python3 analysis.py [-h] [-r TRAJ] [-f TRAJ] 
                             [-s TOPOL] [-n INDEX] [-g GRO] 
                             [-d EDR] [-m MDP] [-p TOP]
                             [-o OUT]

required arguments:
  -r TRAJ              file    .trr trajectory file to analyze 
                                    for lateral pressure profile.
  -f TRAJ              file    .xtc trajectory file to analyze 
                                    for all other properties.
  -s TOPOL             file    .tpr topology file with the same 
                                    atom numbering as the trajectory.
  -n INDEX             file    .ndx index file with three groups:
                                    System, Solvent, Membrane
  -g GRO               file    .gro topology file with the same 
                                    atom numbering as the trajectory.
  -d EDR               file    .edr portable energy file matching 
                                    trajectory.
  -m MDP               file    .mdp molecular dynamics parameters file 
                                    matching trajectory.
  -p TOP               file    .top topology file with the same 
                                    atom numbering as the trajectory.
  -o OUT               dir     output directory for analysis files

optional arguments:
  -h, --help                   display this help message and exit
"""

# Handle command-line arguments
try:
	opts, args = getopt.getopt(sys.argv[1:], "hr:f:s:n:g:d:m:p:o:", 
		["trrFP=", "xtcFP=", "tprFP=", "indFP=", "groFP=", "edrFP=", "mdpFP=", "topFP=", "outDir="])
except getopt.GetoptError as err:
	print(err)
	sys.exit(2)

for opt, arg in opts:
	if opt == '-h':
		print(usage)
		exit(2)

#################################################################################
###############################  START OF SCRIPT  ###############################
#################################################################################

resultFile = open(outDir + 'analysisResults.txt',"w+")
syst = mdreader.MDreader()

# Parse bilayer composition
lipidNameList, leaflet_down, leaflet_up, asym = parse_bilayer_composition()
resultFile.write(f"Asymmetry {asym}\n\n")

# Parse simulation parameters
totalSimulationTime, temperature = parse_md_param()

# Create new index file grouped by lipid headgroup
selString = "(echo del 0-100; echo a PO4 NH3 CNO ROH; echo name 0 headgroups; echo q)"
headgroupsFP = outDir + 'headgroups.ndx'
createIndexFile(selString, indFP, headgroupsFP)

# Separate bilayer by leaflet into two index files using fatslim
leafletsFP = outDir + 'bilayer_leaflet_0000.ndx'
membranes(leafletsFP)
selString1 = "(echo del 0; echo q)"
selString2 = "(echo del 1; echo q)"
downFP = outDir + 'down.ndx'
upFP = outDir + 'up.ndx'
createIndexFile(selString1, leafletsFP, downFP)
createIndexFile(selString2, leafletsFP, upFP)

# Center trajectory files
xtcOut = xtcFP[:-4] + '_centered' + xtcFP[-4:]
trrOut = trrFP[:-4] + '_centered' + trrFP[-4:]
centerTrajFiles(xtcOut, trrOut)
xtcFP, trrFP = xtcOut, trrOut

# Calculate bilayer thickness using fatslim and write to output file
resultFile.write("Bilayer Thickness = %s with a standard error of %s\n" % avg_and_bse(thickness()))

# Calculate area per lipid for each leaflet using fatslim and write to output file
(apl_down, apl_up) = apl()
resultFile.write("Area per Lipid (Lower) = %s with a standard error of %s\n" % avg_and_bse(apl_down))
resultFile.write("Area per Lipid (Upper) = %s with a standard error of %s\n" % avg_and_bse(apl_up))

# Calculate area compressibility for each leaflet and write to output file
areacomp_down = computeAreaComp(apl_down, finalSimulationTime/5, finalSimulationTime, sum(leaflet_down), temperature)
areacomp_up = computeAreaComp(apl_up, finalSimulationTime/5, finalSimulationTime, sum(leaflet_up), temperature)
resultFile.write(f"Area compressibility (Lower) = %s with a standard error of %s\n") % areacomp_down
resultFile.write(f"Area compressibility (Upper) = %s with a standard error of %s\n") % areacomp_up

# Calculate order parameters for each lipid type
for i in range(len(lipidNameList)):
	orderParam_down = []
	orderParam_up = []
	for j in range(1, 5): # using 5 blocks
		startTime = (j / 5) * finalSimulationTime
		endTime = ((j + 1) / 5) * finalSimulationTime
		if leaflet_down[i] > 0:
			down = computeOrderParameter(xtcFP, tprFP, downFP, startTime, endTime, 5, 0, 0, 1, leaflet_down[i], lipidNameList[i])
			orderParam_down.append(down)
		if leaflet_up[i] > 0:
			up = computeOrderParameter(xtcFP, tprFP, upFP, startTime, endTime, 5, 0, 0, 1, leaflet_up[i], lipidNameList[i])
			orderParam_up.append(up)

	# compute avg/sd and write to output
	if leaflet_down[i] > 0:
		ordAvg = sum(orderParam_down) / len(orderParam_down)
		ordErr = statistics.pstdev(orderParam_down) / 2
		resultFile.write(f"Average order parameter (Lower) for {lipidNameList[i]} = {ordAvg} with a standard error of {ordErr}\n")
	if leaflet_up[i] > 0:
		ordAvg = sum(orderParam_up) / len(orderParam_up)
		ordErr = statistics.pstdev(orderParam_up) / 2
		resultFile.write(f"Average order parameter (Upper) for {lipidNameList[i]} = {ordAvg} with a standard error of {ordErr}\n")	

# Collect head groups needed for diffusion index file
groups = []
for i in range(len(lipidNameList)):
	lipid = lipidNameList[i]
	headgroup = get_headgroup(lipid)
	if leaflet_down[i] > 0:
		groups.append(f"1 & r {lipid} & a {headgroup}")
	if leaflet_up[i] > 0:
		groups.append(f"0 & r {lipid} & a {headgroup}")

# Create index file for lateral diffusion coefficients
selString = "(echo name 0 up; echo name 1 down; "		
for g in groups:
	selString += 'echo "' + g + '"; '
selString += "echo q)"
diffusionFP = outDir + 'diffusion.ndx'
createIndexFile(selString, leafletsFP, diffusionFP)

# Calculate diffusion coefficients for each lipid and write to output file
computeDiffusionCoeffs(groups, diffusionFP)

# Calculate cholesterol flip flop rate and average # of CHOL in each leaflet
rate, lowerN, upperN = cholesterolFlipFlop()
resultFile.write(f"CHOL flip flip rate: {rate} \n")
resultFile.write(f"CHOL leaflet counts: Upper {upperN} Lower {lowerN} \n")

# Calculate surface tension and write to output file
top_sf, bot_sf = surface_tension()
resultFile.write(f"Top Leaflet Surface Tension {top_sf} \n")
resultFile.write(f"Bottom Leaflet Surface Tension {bot_sf} \n")

resultFile.flush()
resultFile.close()

#################################################################################
###############################   END OF SCRIPT   ###############################
#################################################################################

"""

Boilerplate functions for property analysis

Some functions use dictionaries defined at the bottom of the file.
These dictionaries store the string representations by lipid type
and can be easily extended to include custom lipids outside the scope of
our current project.

"""

# Return name of head group bead for a lipid
def get_headgroup(lipid):
	lipidID = lipid[2:]
	headgroup = ""
	# Determine headgroup
	if lipidID in ["PC", "PE", "PS", "PG", "PA", "SM", "P6"]:
		headgroup = "PO4"
	elif lipidID == "PI":
		headgroup = "CP"
	elif lipid == "CHOL":
		headgroup = "ROH"
	else:
		print(f"Error in analysis: {lipid} headgroup undefined")
		exit(1)
	return(headgroup)

#given a list of measurements, find their mean and block standard error 
def avg_and_bse(num_data):
	dataMean = sum(num_data)/len(num_data)
	blockLengths = np.linspace(1, len(num_data)/5, num=len(num_data)/5)
	bse = []
	for n in blockLengths:
		numBlocks = len(num_data)/n
		blockAverages = []
		for i in range(int(numBlocks)):
			avgi = sum(num_data[int(i*n):int((i+1)*n)])/n
			blockAverages.append(avgi)
		# now we have means for all our blocks: compute standard error
		if len(blockAverages) > 1:
			blockse = statistics.stdev(blockAverages)/math.sqrt(numBlocks)
		else:
			blockse = 0
		bse.append(blockse)
	bseError = max(bse)
	return(dataMean, bseError)

# generate new index file with given selection string of form
# (SEL1; SEL2; ...) a sequence of selection commands 
def createIndexFile(selString, indexIn, indexOut): 
	command = selString + f" | gmx make_ndx -f {groFP} -n {indexIn} -o {indexOut}"
	subprocess.call(command, shell=True)
	return

# center trajectory files (.xtc, .trr)
# centers system relative to membrane in box
def centerTrajFiles(xtcOut, trrOut):
	# xtc 
	command = f"(echo 2; echo 0) | gmx trjconv -f {xtcFP} -o {xtcOut} -n {indFP} -center -pbc mol -s {tprFP}"	
	subprocess.call(command, shell=True)
	# trr
	command = f"(echo 2; echo 0) | gmx trjconv -f {trrFP} -o {trrOut} -n {indFP} -center -pbc mol -s {tprFP}"	
	subprocess.call(command, shell=True)
	return

# Create index file grouping leaflets separately using fatslim
def membranes(indexOut):
	fp = "/p/lustre1/lamarpat/fatslim/fatslim"
	subprocess.call(["python2.7", fp, "membranes", "-c", groFP,  "-n", headgroupsFP, "--output-index", indexOut])
	return

# Compute bilayer thickness using fatslim 
# Return time series data
def thickness():
	# call fatslim
	fp = "/p/lustre1/lamarpat/fatslim/fatslim"
	xvgOut = outDir + "thickness.xvg"
	subprocess.call(["python2.7", fp, "thickness", "-c", groFP, "-n", headgroupsFP, "-t", xtcFP, "-b", str(finalSimulationTime / 5), "-e", str(finalSimulationTime), "--plot-thickness", xvgOut])
	# parse the output xvg
	metadata, num_data = parse_xvg(xvgOut)
	return(num_data[1])

# Compute area per lipid using fatslim 
# Return time series data for both leaflets (Lower, Upper)
def apl():
	# call fatslim
	fp = "/p/lustre1/lamarpat/fatslim/fatslim"
	xvgOut = outDir + "apl.xvg"
	subprocess.call(["python2.7", fp, "apl", "-c", groFP, "-n", headgroupsFP, "-t", xtcFP, "-b", str(finalSimulationTime / 5), "-e", str(finalSimulationTime), "--plot-apl", xvgOut])
	# parse the output xvg
	metadata, num_data = parse_xvg(outDir + "apl.xvg")
	return(num_data[2], num_data[3])

'''
Compute the projected area per lipid (A_0) and area compressibility modulus (K_A). 
Using "Box-X" and "Box-Y" extracted from an energy file given in the arguments and 
dividing this area with the given number of lipids in each monolayer. 

See for example {Marrink, Vries, Mark, J. Phys. Chem. 2004, 108:750-760} and {Feller
SE and Pastor RW, J. Chem. Phys. 1999, 111:1281-1287}.
Small system will overestimate the K(A)
K(A) = kT<A> / (N<(A - A0)2>)

WARNING: will override/output data to ener-area.xvg, ener-area-A0-block.xvg and ener-area-KA-block.xvg
'''
def computeAreaComp(area_lipid, start, end, lipids_per_mon, temperature):
	
    scaled_k = 0.013806488 # Boltzmann constant (1.3806488 x 10^23 J/K, but here scaled so final K(A) values is in mN/m)
    line_count = len(area_lipid)

    # Get ave +/- sd for area per lipid
    val_A0 = numpy.average(area_lipid)
    val_A0_sd_r = numpy.std(area_lipid)

    # Use block averaging to get the SD and SE of the average area per lipid
    # Use minimum x5 blocks, report SD SE as max value in the last 20% in block curve
    block_max = int(line_count/5)
    block_range = range(1, block_max, 5)
    block_in_top_20 = block_max * 0.8
    val_A0_sd = -1
    val_A0_se = -1
    for block_size in block_range:
        blocks = range(block_size, line_count+1, block_size)
        numBlocks = int(line_count / block_size)
        blocked_A0 = []
        index = 0
        last_block = 0
        for ti in blocks:
            blocked_A0.append(numpy.average(area_lipid[(ti-block_size):ti]))

        currentSTD = numpy.std(blocked_A0)
        currentSE  = currentSTD / math.sqrt(numBlocks)

        if (block_size >=  block_in_top_20):
            if val_A0_sd < currentSTD:
                val_A0_sd = currentSTD
            if val_A0_se < currentSE:
                val_A0_se = currentSE

    # Calculate area compressibility modulus K(A) = kTA0 / N<(A - A0)^2>
    val_ave_area_var = 0
    area_lipid_squear_diff = numpy.zeros(line_count)
    # For each timepoint
    for ti in range(line_count):
        area_lipid_squear_diff[ti] = (area_lipid[ti] - val_A0)**2
        val_ave_area_var += area_lipid_squear_diff[ti]
    # End for
    val_ave_area_var = val_ave_area_var / line_count
    val_scaling = (scaled_k * temperature * val_A0) / lipids_per_mon  # kT A0 / N
    val_KA = val_scaling * (1 / val_ave_area_var)

    # Calculate SD for area compressibility modulus K(A) - need to do block averaging 
    # Use minimum x5 blocks, report SD SE as max value in the last 20% in block curve
    block_max = int(line_count/5)
    block_range = range(1, block_max, 5)
    block_in_top_20 = block_max * 0.8
    val_KA_sd = -1
    val_KA_se = -1
    for block_size in block_range:
        blocks = range(block_size, line_count, block_size)
        numBlocks = int(line_count / block_size)
        blocked_A0 = []
        index = 0
        for ti in blocks:
            current_block = []
            while (index < ti):
                current_block.append(area_lipid_squear_diff[index])
                index += 1
    
            blocked_A0.append(val_scaling * (1 / numpy.average(current_block)))

        currentSTD = numpy.std(blocked_A0)
        currentSE  = currentSTD / math.sqrt(numBlocks)
        print("%7i    %10.8f    %10.8f    %10.8f" % (block_size, numpy.average(blocked_A0), currentSTD, currentSE), file=outFile)

        if (block_size >=  block_in_top_20):
            if val_KA_sd < currentSTD:
                val_KA_sd = currentSTD
            if val_KA_se < currentSE:
                val_KA_se = currentSE

    return(val_KA, val_KA_se)


'''
Compute (second rank) order parameter, defined as:

  P2 = 0.5*(3*<cos²(theta)> - 1)

where "theta" is the angle between the bond and the bilayer normal.
P2 = 1      perfect alignement with the bilayer normal
P2 = -0.5   anti-alignement
P2 = 0      random orientation

WARNING function will output all frames in one go, into files called frame_dump_XXX.gro and 
then remove them so may cause collisions with other processes.
'''
def computeOrderParameter(trajfile, tprfile, indexfile, initial_time, final_time,
	                      traj_skip, normalx, normaly, normalz, number_of_lipids, lipid_type):
	# (normalized) orientation of bilayer normal
	orientation_of_bilayer_normal = [normalx, normaly, normalz]
	norm = sqrt(orientation_of_bilayer_normal[0]**2 + orientation_of_bilayer_normal[1]**2 + orientation_of_bilayer_normal[2]**2)
	for i in range(3):
		orientation_of_bilayer_normal[i] /= norm
		stdout.write("(Normalized) orientation of bilayer normal: ( %.3f | %.3f | %.3f ).\n" % (
			orientation_of_bilayer_normal[0], \
			orientation_of_bilayer_normal[1], \
			orientation_of_bilayer_normal[2]  \
			))

	# lookup bonds for lipid type
	bond_names = bondDict[lipid_type]

	# Output all frame using trjconv 
	command = "echo %s | gmx trjconv -f %s -s %s -n %s -b %i -e %i -sep -skip %i -pbc whole -o %sframe_dump_.gro > /dev/null" % (lipid_type, trajfile, tprfile, indexfile, initial_time, final_time, traj_skip, outDir)
	subprocess.call(command, shell=True)

	# For each dumped frame
	order_parameters = []
	file_count = 0
	bonds = []
	while True:
		filename = outDir + "frame_dump_" + str(file_count) + ".gro"
		if not path.isfile(filename) or path.getsize(filename) == 0:
			break

		# compute order parameter for each bond, for each snapshot
		current_order_parameters = []
		# bonds respectively involved in the head,
		#                             in the junction head-tail,
		#                             in each tail
		bonds = []

		for bond_name in bond_names.split():
			bonds.append(bond_name.split("-"))

		for bond in bonds:

			# parse .gro file, grep bead coordinates
			first, second = read_gro(filename, bond)

			# compute order parameter for each lipid
			order_parameter = 0.0
			for i in range(number_of_lipids):
				# vector between the two previous beads (orientation doesn't matter)
				vector = [0.0, 0.0, 0.0]
				for j in range(3):
					vector[j] = first[i][j] - second[i][j]
					norm2 = vector[0]**2 + vector[1]**2 + vector[2]**2
				# compute projection on the bilayer normal
				projection = vector[0]*orientation_of_bilayer_normal[0] + vector[1]*orientation_of_bilayer_normal[1] + vector[2]*orientation_of_bilayer_normal[2]
				# order parameter
				order_parameter += projection**2/norm2

			# compute final averaged order parameter
			# store everything in lists
			current_order_parameters.append(0.5*(3.0*(order_parameter/number_of_lipids) - 1.0))
		order_parameters.append(current_order_parameters)


		remove(filename)
		file_count += 1
	# End of while loop

	# average order parameter
	averaged_order_parameters = []
	for i in range(len(bonds)):
		sum = 0.0
		for j in range(len(order_parameters)):
			sum += order_parameters[j][i]
			averaged_order_parameters.append(sum/len(order_parameters))

	# Calculate abs average order parameters <Sn> (for carbon chains only)
	ave_chain_s = 0
	for i in averaged_order_parameters[3:]: 
		ave_chain_s += abs(i)

	return(ave_chain_s / (len(averaged_order_parameters)-3))

# Calculate the diffusion coefficients for all groups in the index file
# using gromacs and write results to output file
def computeDiffusionCoeffs(groups, indexFP):
	selString = "("
	for i in range(len(groups)):
		selString += f"echo {i+2}; " # first two groups in index file should be the two leaflets
	selString = selString[:-2] + ")"
	command = f"gmx msd -s {tprFP} -f {xtcFP} -lateral z -n {indexFP} -ngroup {len(groups)} -o {outDir}diffusion_coeffs.xvg -b {finalSimulationTime / 5} -e {finalSimulationTime}"
	subprocess.call(selString + " | " + command, shell=True)

	# parse output xvg files for coefficients 
	diffusionFile = open(f"{outDir}diffusion_coeffs.xvg", "r")
	for line in diffusionFile.readlines():
		if line[0:4] == "# D[":
			tokens = line[4:].split("_&_")
			layer = tokens[0]
			lipid = tokens[0]
			coeff = line.split("=")[1].split("(")[0]
			err = line.split("+/-")[1].split(")")[0]
			if layer == "up":
				resultFile.write("Lateral Diffusion Coef (Upper) for %s = %s with a standard error of %s\n" % (lipid, coeff, err))
			else:
				resultFile.write("Lateral Diffusion Coef (Lower) for %s = %s with a standard error of %s\n" % (lipid, coeff, err))
	return

# Simple gaussian function for curve fitting
def gaussian(x, mean, sd, a):
    norm = a * numpy.exp(-(x - mean)**2 / (2 * sd**2))
    return(norm)

# Fit x2 gaussians
def doublegfit(p, y, x):
    m1, m2, sd1, sd2, a1, a2 = p
    if m1 < 0 or m2 < 0 or m1 > 1000 or m2 > 1000 or a1 > 5 or a2 > 5 or a1 < 0 or a2 < 0:
        return(1e10)
    else:
        y_fit = gaussian(x, m1, sd1, a1) + gaussian(x, m2, sd2, a2)
        return(y - y_fit)

def get_rohs_middle():
    global plsq
    cArray = allPO4.positions[:,2]
    y_raw,x_temp = numpy.histogram(cArray, bins=binRlenBinCounts, range=(0,binRlenMax)) # 0 is data, 1 is x spacing
    y = y_raw.astype(numpy.float32)/numpy.sum(y_raw)
    plsq = leastsq(doublegfit, plsq, args = (y,x))[0]
    return(rohs.positions[:,2] - (plsq[1] + plsq[0]) / 2)

# Calculate the average number of cholesterol molecules in each leaflet and the average flip flop rate
def cholesterolFlipFlop():
	# Bead has to be above or bellow this distance from center to be in on or the other bilayers
	fliplim = 8

	rohs = syst.select_atoms("name ROH")
	allPO4 = syst.select_atoms("name PO4")

	binRlenBinCounts = 800   
	binRlenMax = 200 

	# Initialize distribution parameters using statistics over full simulation
	i = int(len(allPO4.positions) / 2)
	top = allPO4.positions[:i,2]
	bot = allPO4.positions[i:,2]
	topZ = float(sum(top) / len(top))
	botZ = float(sum(bot) / len(bot))
	topSD = float(math.sqrt(sum([pow(t - topZ, 2) for t in top]) / (len(top) - 1)))
	botSD = float(math.sqrt(sum([pow(t - botZ, 2) for t in bot]) / (len(bot) - 1)))

	# Initial guesses for x2 gaussian used in leastsq 
	# (mean-inner, mean-outer, sd-inner, sd-outer, amplitude-inner, ampitude-outer) 
	plsq = numpy.array([botZ, topZ, botSD, topSD, 0.1, 0.1])    

	x = numpy.linspace(0,binRlenMax,binRlenBinCounts)

	cdx = numpy.array(syst.do_in_parallel(get_rohs_middle))

	thresh = (cdx>fliplim).astype(numpy.int)
	thresh -= cdx<-fliplim
	thresh = thresh.T

	thresh = numpy.column_stack((thresh, numpy.ones(len(thresh))*1000))
	thresh = thresh[thresh!=0] # Becomes 1D
	diff = numpy.diff(thresh)

	upflips = (diff==2).sum()
	downflips = (diff==-2).sum()
	numflip = upflips + downflips


	# Also count number in each leaflet and report 
	thresh = (cdx>fliplim).astype(numpy.int)
	countUpper = numpy.sum(thresh).astype(numpy.float)/syst.totalframes
	thresh = (cdx<-fliplim).astype(numpy.int)
	countLower = numpy.sum(thresh).astype(numpy.float)/syst.totalframes

	return(numflip / syst.totalframes, countLower, countUpper)

# Calculate Lateral Pressure Profile using modified gromacs
def surface_tension():
	# Create a new .tpr file in the modified version of gromacs
	command = f'gmx_LS grompp -f {mdoFP} -c {groFP} -n {indFP} -p {topFP} -o {runDir}/20fs/20fs_LS.tpr'
	subprocess.call(command, shell=True)

	# Translate the bilayer before centering
	#command = f'gmx editconf -f {trrFP} -o {runDir}/20fs/20fs_translated.trr -translate 0 0 2.5'
	#subprocess.call(command, shell=True)

	# Center the bilayer in the simulation box
	command = f'(echo 2; echo 0) | gmx_LS trjconv -f {trrFP} -o {runDir}/20fs/20fs_centered.trr -n {indFP} -center -pbc mol -s {runDir}/20fs/20fs_LS.tpr -b 8000000'
	subprocess.call(command, shell=True)

	# Create a new directory to hold all of the stress files
	command = f"mkdir {runDir}/20fs/stress"
	subprocess.call(command, shell=True)
	command = f"cd {runDir}/20fs/stress"
	subprocess.call(command, shell=True)

	# Rerun the trajectory file to calculate the 3-d presure profile
	command = f"gmx_LS mdrun -s {runDir}/20fs/20fs_LS.tpr -rerun {runDir}/20fs/20fs_centered.trr -g {runDir}/20fs/stress/md.log -e {runDir}/20fs/stress/ener.edr -o {runDir}/20fs/stress/traj.trr -ols {runDir}/20fs/stress/localstress -localsgrid 0.3 -lsfd gld"
	subprocess.call(command, shell=True)

	# Calculate the Pressure Tensor and save it as localstress.dat0
	command =f"tensortools -f {runDir}/20fs/stress/localstress.dat0 --prof z -o {runDir}/20fs/stress/stress_z.txt"                
	subprocess.call(command, shell=True)

	# Initiate the dictionary for P_n, z, xx , yy , where P_N = -(zz)
	pressure_norm, stress = [], []
	x_x, y_y, zval = [], [], []

	# Read in stress tensor and calculate lateral pressure profile
	stress = np.loadtxt(strFP + '/stress_z.txt')
	a, x, y, z = [], [], [], []
	for j in range(len(stress)):
		a = np.append(a,-1*stress[j][9])
		x = np.append(x, stress[j][1])
		y = np.append(y, stress[j][5])
		z = np.append(z, stress[j][0])
        
	pressure_norm, x_x, y_y, z_val = [a], [x], [y], [z]
   	
   	# Calculating P_L where P_L = -(xx+yy)/2
	pressure_lat = -1*(x_x[0] + y_y[0]) /2

	# Calculate Lateral Pressure = P_L-P_N
	lat_pres_prof = pressure_lat-pressure_norm
	
	# Now spline the data in order to calculate the lateral pressure profile
	z_val_1 = z_val[0][0:int(len(z_val[0])/2)]
	z_val_2 = z_val[0][int(len(z_val[0])/2):]

	lat_pres_prof_1 = lat_pres_prof[0][0:int(len(lat_pres_prof[0])/2)]
	lat_pres_prof_2 = lat_pres_prof[0][int(len(lat_pres_prof[0])/2):]

	s_bot = UnivariateSpline(z_val_1, lat_pres_prof_1, s=1)
	xs_bot = linspace(z_val[0][0], z_val[0][int(len(z_val[0])/2)], 50 )
	ys_bot = s_bot(xs_bot)

	s_top = UnivariateSpline(z_val_2, lat_pres_prof_2, s=1)
	xs_top = linspace(z_val[0][int(len(z_val[0])/2)], z_val[0][-1], 50 )
	ys_top = s_top(xs_top)

	# Now we calculate the leaflet surface tension
	top_sf = np.trapz(ys_bot)
	bot_sf = np.trapz(ys_top)

	return(top_sf, bot_sf)

"""

Auxillary functions for parsing and reading from 
gromacs and analysis files

"""

# Parse .mdp file for simulation parameters
# Return total simulation time
def parse_md_param():
	# read values out of mdp file
	mdpFile = open(mdpFP, "r")
	precision = 0
	dt = 0
	nsteps = 0
	temperature = 0
	for line in mdpFile.readlines():
		if line[0:2] == "dt":
			dt = float(line.split("=")[1])
		if line[0:6] == "nsteps":
			nsteps = int(line.split("=")[1])
		if line[0:22] == "compressed-x-precision":
			precision = int(line.split("=")[1])
		if line[0:8] == "gen_temp":
			temperature = int(line.split("=")[1])
	timeBetweenFrames = precision * dt 
	finalSimulationTime = dt * nsteps
	mdpFile.close()

	return finalSimulationTime, temperature

# Parse .top file generated by INSANE 
# Return bilayer composition and asymmetry
def parse_bilayer_composition():
	# read values out of topology file
	lipidNameList = []
	leaflet1 = []
	leaflet2 = []
	topFile = open(topFP, "r")
	p = re.compile('[A-Z]{3}[A-Z\d] ')
	for line in topFile.readlines():
		if bool(p.match(line)):
			if line[0:4] in lipidNameList:
				i = lipidNameList.index(line[0:4])
				if len(leaflet2) < i:
					while len(leaflet2) < i:
						leaflet2.append(0)
				leaflet2.append(int(line[4:]))
			else:
				lipidNameList.append(line[0:4])
				leaflet1.append(int(line[4:]))
	while len(leaflet2) < len(leaflet1):
		leaflet2.append(0)
	asym = abs(sum(leaflet2) - sum(leaflet1))
	topFile.close()

	if sum(leaflet2) > sum(leaflet1):
		upper = leaflet2
		lower = leaflet1
	else:
		upper = leaflet1
		lower = leaflet2

	return (lipidNameList, lower, upper, asym)

# Functions for handling xvg files
def parse_xvg(fname, sel_columns='all'):
    """Parses XVG file legends and data"""
    _ignored = set(('legend', 'view'))
    _re_series = re.compile('s[0-9]+$')
    _re_xyaxis = re.compile('[xy]axis$')

    metadata = {}
    num_data = []
    
    metadata['labels'] = {}
    metadata['labels']['series'] = []

    ff_path = os.path.abspath(fname)
    if not os.path.isfile(ff_path):
        raise IOError('File not readable: {0}'.format(ff_path))
    
    with open(ff_path, 'r') as fhandle:
        for line in fhandle:
            line = line.strip()
            if line.startswith('@'):
                tokens = shlex.split(line[1:])
                if tokens[0] in _ignored:
                    continue
                elif tokens[0] == 'TYPE':
                    if tokens[1] != 'xy':
                        raise ValueError('Chart type unsupported: \'{0}\'. Must be \'xy\''.format(tokens[1]))
                elif _re_series.match(tokens[0]):
                    metadata['labels']['series'].append(tokens[-1])
                elif _re_xyaxis.match(tokens[0]):
                    metadata['labels'][tokens[0]] = tokens[-1]
                elif len(tokens) == 2:
                    metadata[tokens[0]] = tokens[1]
                else:
                    print('Unsupported entry: {0} - ignoring'.format(tokens[0]), file=sys.stderr)
            elif line[0].isdigit():
                num_data.append(map(float, line.split()))
    
    num_data = list(zip(*num_data))

    if not metadata['labels']['series']:
        for series in range(len(num_data) - 1):
            metadata['labels']['series'].append('')

    # Column selection if asked
    if sel_columns != 'all':
        sel_columns = map(int, sel_columns)
        x_axis = num_data[0]
        num_data = [x_axis] + [num_data[col] for col in sel_columns]
        metadata['labels']['series'] = [metadata['labels']['series'][col - 1] for col in sel_columns]
    
    return metadata, num_data

# parse a .gro file
# return a list of coordinates
def read_gro(file, atoms):
	line_counter = 0
	number_of_particles = 0
	first, second = [], []
	for line in open(file):
		if line_counter == 1:
			number_of_particles = int(line)
		elif line_counter > 1 and line_counter < number_of_particles + 2:
			if line[10:15].strip() == atoms[0]:
				first.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
			elif line[10:15].strip() == atoms[1]:
				second.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
		line_counter += 1
	return [first, second]

"""

Lipid Dictionaries

"""

phosphatidylcholine_bond_names = " NC3-PO4 PO4-GL1 GL1-GL2 "
phosphatidylethanolamine_bond_names = " NH3-PO4 PO4-GL1 GL1-GL2 "
phosphatidylserine_bond_names = " CNO-PO4 PO4-GL1 GL1-GL2 "
phosphatidylglycerol_bond_names = " GL0-PO4 PO4-GL1 GL1-GL2 "
phosphatidic_acid_bond_names = " PO4-GL1 GL1-GL2 "
phosphatidylinositol_bond_names = " C1-C2 C1-C3 C2-C3 C1-PO4 PO4-GL1 GL1-GL2 "

bondDict = {
	# PCs
	"DAPC":phosphatidylcholine_bond_names + "GL1-D1A GL2-D1B D1A-D2A D2A-D3A D3A-D4A D4A-C5A D1B-D2B D2B-D3B D3B-D4B D4B-C5B",
	"DBPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B",
	"DFPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A GL2-C1B C1B-D2B D2B-D3B D3B-D4B",
	"DGPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B",
	"DIPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-D2B D2B-D3B D3B-C4B",
	"DLPC":phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B",
	"DNPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B",
	"DOPC":phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-D2A D2A-C3A C3A-C4A C1B-D2B D2B-C3B C3B-C4B",
	"DPPC":phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B",
	"DRPC":phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-D6B",
	"DTPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A GL2-C1B C1B-C2B",
	"DVPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-D3B D3B-C4B",
	"DXPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B",
	"DYPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-C3A GL2-C1B C1B-D2B D2B-C3B",
	"LPPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B",
	"PAPC":phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PEPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PGPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PIPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"POPC":phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PRPC":phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PUPC":phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	# PEs
	"DAPE":phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-C5B",
	"DBPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B",
	"DFPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A GL2-C1B C1B-D2B D2B-D3B D3B-D4B",
	"DGPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B",
	"DIPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-D2B D2B-D3B D3B-C4B",
	"DLPE":phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B",
	"DNPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B",
	"DOPE":phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-D2A D2A-C3A C3A-C4A C1B-D2B D2B-C3B C3B-C4B",
	"DPPE":phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B",
	"DRPE":phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-D6B",
	"DTPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A GL2-C1B C1B-C2B",
	"DUPE":phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B",
	"DVPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-D3B D3B-C4B",
	"DXPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B",
	"DYPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-C3A GL2-C1B C1B-D2B D2B-C3B",
	"LPPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B",
	"PAPE":phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PGPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PIPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"POPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PQPE":phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PRPE":phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PUPE":phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	# PSs
	"DAPS":phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-C5B",
	"DBPS":phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B",
	"DFPS":phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A GL2-C1B C1B-D2B D2B-D3B D3B-D4B",
	"DGPS":phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B",
	"DIPS":phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-D2B D2B-D3B D3B-C4B",
	"DLPS":phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A GL2-C1B C1B-C2B C2B-C3B",
	"DNPS":phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B",
	"DOPS":phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-D2B D2B-C3B C3B-C4B",
	"DPPS":phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"DRPS":phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-D6B",
	"DTPS":phosphatidylserine_bond_names + "GL1-C1A C1A-C2A GL2-C1B C1B-C2B",
	"DUPS":phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B",
	"DVPS":phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-D3B D3B-C4B",
	"DXPS":phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B",
	"DYPS":phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-C3A GL2-C1B C1B-D2B D2B-C3B",
	"LPPS":phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B",
	"PAPS":phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PGPS":phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PIPS":phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"POPS":phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PQPS":phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PRPS":phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PUPS":phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	# PGs
	"DAPG": phosphatidylglycerol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-C5B",
	"DBPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B",
	"DFPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A GL2-C1B C1B-D2B D2B-D3B D3B-D4B",
	"DGPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B",
	"DIPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-D2B D2B-D3B D3B-C4B",
	"DLPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A GL2-C1B C1B-C2B C2B-C3B",
	"DNPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B",
	"DOPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-D2B D2B-C3B C3B-C4B",
	"DPPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"DRPG": phosphatidylglycerol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-D6B",
	"DTPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A GL2-C1B C1B-C2B",
	"DVPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-D3B D3B-C4B",
	"DXPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B",
	"DYPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-C3A GL2-C1B C1B-D2B D2B-C3B",
	"LPPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B",
	"PAPG": phosphatidylglycerol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PGPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PIPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"POPG": phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PRPG": phosphatidylglycerol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	# PAs
	"DAPA": phosphatidic_acid_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-C5B",
	"DBPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B",
	"DFPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A GL2-C1B C1B-D2B D2B-D3B D3B-D4B",
	"DGPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B",
	"DIPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-D2B D2B-D3B D3B-C4B",
	"DLPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A GL2-C1B C1B-C2B C2B-C3B",
	"DNPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B",
	"DOPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-D2B D2B-C3B C3B-C4B",
	"DPPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"DRPA": phosphatidic_acid_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-D6B",
	"DTPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A GL2-C1B C1B-C2B",
	"DVPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-D3B D3B-C4B",
	"DXPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B",
	"DYPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-C3A GL2-C1B C1B-D2B D2B-C3B",
	"LPPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B",
	"PAPA": phosphatidic_acid_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PGPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PIPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"POPA": phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PRPA": phosphatidic_acid_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PUPA": phosphatidic_acid_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	# PIs
	"DPPI": phosphatidylinositol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PAPI": phosphatidylinositol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PIPI": phosphatidylinositol_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"POPI": phosphatidylinositol_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"PUPI": phosphatidylinositol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B",
	# Other Lipids
	"DPSM": "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A T1A-C2A C2A-C3A AM2-C1B C1B-C2B C2B-C3B C3B-C4B",
	"CHOL": "C1-C2 ROH-R2 ROH-R3 R2-R3 C1-R2 C1-R3",
	"PAP6": "C1-C2 C1-C3 C2-C3 C3-P4 C3-P5 C2-P5 C2-P4 C1-PO4 PO4-GL1 GL1-GL2 GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B"
}

