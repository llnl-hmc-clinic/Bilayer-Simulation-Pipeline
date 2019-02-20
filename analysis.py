import numpy as np 
import math
import os
import do_order_module
#import property_plot
import area_lipid_g5
import subprocess
import re
import statistics
import shlex

'''
Structure for property lists for area per lipid, area compressibility, and 
bilayer thickness is as follows:
    propList = [asymmetry1, asymmetry2, asymmetry3, ...] asymmetries in ascending order
    where each asymmetry holds the average calculated property value at that asymmetry

Structure for property lists for order parameter and lateral diffusion is as follows:
    propList = [asymmetry1, asymmetry2, asymmetry3, ...] asymmetries in ascending order
    asymmetry = [lipid1, lipid2, ...] each asymmetry contains a list of values for each lipid
    where each lipid holds the average calculated property value for that lipid at that asym

 The structure for the error lists is similar to their respective property lists, except 
 holding standard errors instead of average values

'''
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
    print(num_data)

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

def avg_and_bse(num_data):
	''' given a list of measurements, find their mean and find their block standard
	    error '''
	dataMean = sum(num_data)/len(num_data)
	blockLengths = np.linspace(1, len(num_data), num=len(num_data))
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
	return dataMean, bseError



# input a list of lipid names in the same order as they appear in description.txt, 
# the total time of the simulation (in ps) and the time skip between frames that are saved
# (again in ps)
def analyze():
	# Start from top level directory for a set of runs at different asymmetries
	# Find each individual run
	directoryList = os.listdir()
	filteredList = list(filter((lambda x : x[0:3] == "run"), directoryList))

	# Make a directory in which to put our results, and open a text file in it
	os.mkdir("analysis")
	resultFile = open("analysis/analysisResults.txt","w+")

	# for each run, find the asymmetry value, the 20 fs tpr, gro, xtc, edr, and index file
	for runDir in filteredList:
		groFP = runDir + "/20fs.gro"
		tprFP = runDir + "/20fs.tpr"
		xtcFP = runDir + "/20fs.xtc"
		edrFP = runDir + "/20fs.edr"
		indFP = runDir + "/index.ndx"
		topFP = runDir + "/topol.top"
		
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
		asym = abs(sum(leaflet2) - sum(leaflet1))
		resultFile.write(f"Asymmetry {asym}\n\n")
		topFile.close()

		# read values out of mdp file
		mdpFile = open(runDir + "/mdout.mdp", "r")
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

		# now calculate property values
		'''
		# find order parameter for each lipid type
		for i in range(len(lipidNameList)):
			orderParam = []
			for j in range(1, 5):
				orderParam.append(do_order_module.computeOrderParameter(xtcFP, tprFP, (j / 5) * finalSimulationTime, ((j + 1) / 5) * finalSimulationTime, 5, 0, 0, 1, leaflet1[i] + leaflet2[i], lipidNameList[i]))	
			ordAvg = sum(orderParam) / len(orderParam)
			ordErr = statistics.pstdev(orderParam) / 2
			resultFile.write(f"Average order parameter for {lipidNameList[i]} = {ordAvg}\n")
			resultFile.write(f"Standard error for {lipidNameList[i]} = {ordErr}\n")
		resultFile.write("\n")		

		'''
		# find area compressibility
		# This script gives a standard error as well as a value, so we don't need to break up
		# the simulation by time
		areaComp, areaCompSE = area_lipid_g5.calculate_area_compressibility(edrFP, runDir, finalSimulationTime/5, finalSimulationTime, sum(leaflet2), temperature)
		resultFile.write(f"Area compressibility = {areaComp} with a standard error of {areaCompSE}\n")		
		'''
		# find lateral diffusion with gmx msd. Requires knowing head groups for lipids, which
		# is calculated below
		headGroups = []
		for lipid in lipidNameList:
			if lipid[2:4] == "PC" or lipid[2:4] == "PE" or lipid[2:4] == "PS" or lipid[2:4] == "PG" or lipid[2:4] == "PA" or lipid[2:4] == "SM" or lipid[2:4] == "P6":
			    headGroups.append(f"r {lipid} & a PO4")
			elif lipid[2:4] == "PI":
				headGroups.append(f"r {lipid} & a CP")
			elif lipid[2:4] == "P1":
				headGroups.append(f"r {lipid} & a P1")
			elif lipid[2:4] == "P2":
				headGroups.append(f"r {lipid} & a P2")
			elif lipid[2:4] == "P3":
				headGroups.append(f"r {lipid} & a P3")
			elif lipid == "CHOL":
				headGroups.append(f"r {lipid} & a ROH")
			else:
				print(f"error: {lipid} not of a type included in this script")
		print(lipidNameList)
		#headGroups.append('q')
		print(headGroups)
		
		command = "("		
		for group in headGroups:
			command += 'echo "' + group + '"; '
		command += "echo q) | " + f"gmx make_ndx -f {tprFP} -o analysis/lipid_headgroups.ndx"
		subprocess.call(command, shell=True)

		command = "("
		for i in range(len(headGroups)):
			command += f"echo {i+12}; "
		command = command[:-2] + f") | gmx msd -s {tprFP} -f {xtcFP} -lateral z -n analysis/lipid_headgroups.ndx -ngroup {len(headGroups)} -o analysis/diffusion_coeffs.xvg -b {finalSimulationTime / 5} -e {finalSimulationTime}"
		print(command)		
		subprocess.call(command, shell=True)


		# find diffusion values from the output xvg files and write them to our output file
		metadata, num_data = parse_xvg("analysis/diffusion_coeffs.xvg")
		msdmeta, msddata = running_average(num_data, metadata)
                # TODO: calculate error
		resultFile.write("Lateral Diffusion Coef = %s with a standard error of %s\n" % (msddata[-1], 42))
		'''
		# find area per lipid and bilayer thickness		
		# create index file
		indFP = "analysis/index.ndx"
		command = "(echo del 0-100; echo a PO4 NH3 CNO ROH; echo name 0 headgroups; echo q) | gmx make_ndx -f " + groFP + " -o " + indFP
		subprocess.call(command, shell=True)
			

		
		# call FatSlim commands
		#subprocess.call(["fatslim", "membranes", "-c", groFP,  "-n", indFP, "--output-index", "analysis/bilayer_leaflet.ndx"])
		subprocess.call(["fatslim", "thickness", "-c", groFP, "-n", indFP, "-t", xtcFP, "-b", str(finalSimulationTime / 5), "-e", str(1000 + finalSimulationTime / 5), "--plot-thickness", "analysis/thickness.xvg"])
		subprocess.call(["fatslim", "apl", "-c", groFP, "-n", indFP, "-t", xtcFP, "-b", str(finalSimulationTime / 5), "-e", str(1000 + finalSimulationTime / 5), "--plot-apl", "analysis/apl.xvg", "--plot-area", "analysis/area.xvg"])

		# parse the output xvg
		metadata, num_data = parse_xvg("analysis/thickness.xvg")
		thickness, thicknesserror = avg_and_bse(num_data[1])
		resultFile.write("Bilayer Thickness = %s with a standard error of %s\n" % (thickness, thicknesserror))

		metadata, num_data = parse_xvg("analysis/apl.xvg")
		apl, aplerror = avg_and_bse(num_data[1])
		resultFile.write("Area per Lipid = %s with a standard error of %s\n" % (apl, aplerror))	
			
		# call cholesterol flip flop script
		numUp = 0
		numDown = 0
		for i in range(1, 4):
			subprocess.call(["python3", "trxcounter-bilayer-center_p2.py", "-f", xtcFP, "-s", tprFP, "-o", "data.xvg", "-b", str(i * finalSimulationTime/10), "-e", str((i+1) * finalSimulationTime/10)])
			cholFile = open("dist-chol-all.out", "r")
			line_counter = 0
			for line in cholFile.readlines():
				line_counter += 1
				if line_counter == 2:
					data = line.split("(")[1].split("down)")[0].split(" up/ ")
					numUp += int(data[0])
					numDown += int(data[1])
			cholFile.close()
		resultFile.write("%s flip-flops per ns (%s up/ %s down) \n" % ((numUp + numDown) / 3000, numUp, numDown))

		'''

		# call lateral pressure profile
		#echo 0 | trjconv__LS -f trrFP -o traj_centered.trr -n indFP -center -s tprFP
		#mdrun_LS -s tprFP -rerun traj_centered.trr
		#tensortools -f file.dat0 --prof z -o analysis/stress.txt
		

		'''


