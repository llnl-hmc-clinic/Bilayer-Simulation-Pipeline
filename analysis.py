import numpy as np 
import math
import os, sys
import do_order_module
#import property_plot
import area_lipid_g5
import subprocess
import re
import statistics
import shlex
import multiprocessing

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

def running_average(data, metadata, window=10):
    """
    Performs a running average calculation over all series in data.
    Assumes the first series is the x-axis.
    Appends the series and a new label to the original data and label arrays.
    """

    weights = np.repeat(1.0, window)/window
    s_labels = metadata['labels']['series']
    for n_series, series in enumerate(data[1:]):
        series_rav = np.convolve(series, weights, 'valid')
        s_labels.append('{0} (Av)'.format(s_labels[n_series]))
        data.append(series_rav)
    return metadata, data


def avg_and_bse(num_data):
	''' given a list of measurements, find their mean and find their block standard
	    error '''
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
	return dataMean, bseError



# input a list of lipid names in the same order as they appear in description.txt, 
# the total time of the simulation (in ps) and the time skip between frames that are saved
# (again in ps)
def analyze(dirList, resultsFP):
	numfailures = 0
	resultFile = open(resultsFP,"w+")

	# for each run, find the asymmetry value, the 20 fs tpr, gro, xtc, edr, mdout, trr  and index file
	for runDir in dirList:
		groFP = runDir + "/20fs/20fs.gro"
		tprFP = runDir + "/20fs/20fs.tpr"
		xtcFP = runDir + "/20fs/20fs.xtc"
		edrFP = runDir + "/20fs/20fs.edr"
		indFP = runDir + "/20fs/index.ndx"
		topFP = runDir + "/20fs/topol.top"
		mdoFP = runDir + "/20fs/mdout.mdp"
		trrFP = runDir + "/20fs/20fs.trr"

                
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
		resultFile.write(f"Asymmetry {asym}\n\n")
		topFile.close()

		# read values out of mdp file
		mdpFile = open(runDir + "/20fs/mdout.mdp", "r")
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

		# find area per lipid and bilayer thickness		
		# create index file
		subprocess.call(["mkdir", runDir + "/analysis"])
		headgroupsFP = runDir + "/analysis/headgroups.ndx"
		command = "(echo del 0-100; echo a PO4 NH3 CNO ROH; echo name 0 headgroups; echo q) | gmx make_ndx -f " + groFP + " -o " + headgroupsFP
		subprocess.call(command, shell=True)	
		
		# reduce size of xtc 
		subprocess.call(["gmx", "trjconv", "-f", xtcFP, "-o", runDir + '/20fs/20fs_reduced.xtc', "-trunc", str(4 * finalSimulationTime / 5)])

		xtcFP = runDir + '/20fs/20fs_reduced.xtc'

		command = '(echo 2; echo 0) | gmx trjconv -f ' + xtcFP + ' -o ' + runDir + '/20fs/20fs-center.xtc -center -pbc mol -s ' + tprFP + ' -n ' + indFP 
		subprocess.call(command, shell=True)

		xtcFP = runDir + '/20fs/20fs-center.xtc'	
		
		# call FatSlim commands
		fp = "/p/lustre1/lamarpat/fatslim/fatslim"
		subprocess.call(["python2.7", fp, "membranes", "-c", groFP,  "-n", headgroupsFP, "--output-index", runDir + "/analysis/bilayer_leaflet.ndx"])	
		subprocess.call(["python2.7", fp, "thickness", "-c", groFP, "-n", headgroupsFP, "-t", xtcFP, "-b", str(finalSimulationTime / 5), "-e", str(finalSimulationTime), "--plot-thickness", runDir + "/analysis/thickness.xvg", "--nthreads", "12"])
		subprocess.call(["python2.7", fp, "apl", "-c", groFP, "-n", headgroupsFP, "-t", xtcFP, "-b", str(finalSimulationTime / 5), "-e", str(finalSimulationTime), "--plot-apl", runDir + "/analysis/apl.xvg", "--plot-area", runDir + "/analysis/area.xvg", "--nthreads", "12"])
		
		# parse the output xvg
		metadata, num_data = parse_xvg(runDir + "/analysis/thickness.xvg")
		thickness, thicknesserror = avg_and_bse(num_data[1])
		resultFile.write("Bilayer Thickness = %s with a standard error of %s\n" % (thickness, thicknesserror))

		metadata, num_data = parse_xvg(runDir + "/analysis/apl.xvg")
		apl_down, aplerror_down = avg_and_bse(num_data[2])
		apl_up, aplerror_up = avg_and_bse(num_data[3])
		resultFile.write("Area per Lipid (Lower) = %s with a standard error of %s\n" % (apl_down, aplerror_down))	
		resultFile.write("Area per Lipid (Upper) = %s with a standard error of %s\n" % (apl_up, aplerror_up))
		

		# find area compressibility
		areaComp_down, areaCompSE_down = area_lipid_g5.calculate_area_compressibility(num_data[2], runDir + "/analysis", finalSimulationTime/5, finalSimulationTime, sum(leaflet1), temperature)
		areaComp_up, areaCompSE_up = area_lipid_g5.calculate_area_compressibility(num_data[3], runDir + "/analysis", finalSimulationTime/5, finalSimulationTime, sum(leaflet2), temperature)
		resultFile.write(f"Area compressibility (Lower) = {areaComp_down} with a standard error of {areaCompSE_down}\n")
		resultFile.write(f"Area compressibility (Upper) = {areaComp_up} with a standard error of {areaCompSE_up}\n")


		# create index files for order param
		command = f"(echo del 1; echo q) | gmx make_ndx -f {tprFP} -n {runDir}/analysis/bilayer_leaflet_0000.ndx -o {runDir}/analysis/up.ndx"
		subprocess.call(command, shell=True)
		command = f"(echo del 0; echo q) | gmx make_ndx -f {tprFP} -n {runDir}/analysis/bilayer_leaflet_0000.ndx -o {runDir}/analysis/down.ndx"
		subprocess.call(command, shell=True)

		# find order parameter for each lipid type
		for i in range(len(lipidNameList)):
			if leaflet2[i] > 0:
				try:
					orderParam = []
					for j in range(1, 5):
						up = do_order_module.computeOrderParameter(runDir + "/analysis/", xtcFP, tprFP, runDir + "/analysis/up.ndx", (j / 5) * finalSimulationTime, ((j + 1) / 5) * finalSimulationTime, 5, 0, 0, 1, leaflet2[i], lipidNameList[i])
						orderParam.append(up)
					ordAvg = sum(orderParam) / len(orderParam)
					ordErr = statistics.pstdev(orderParam) / 2
					resultFile.write(f"Average order parameter (Upper) for {lipidNameList[i]} = {ordAvg} with a standard error of {ordErr}\n")
				except:
					numfailures += 1
			if leaflet1[i] > 0:
				try:
					orderParam = []
					for j in range(1, 5):
						down = do_order_module.computeOrderParameter(runDir + "/analysis/", xtcFP, tprFP, runDir + "/analysis/down.ndx", (j / 5) * finalSimulationTime, ((j + 1) / 5) * finalSimulationTime, 5, 0, 0, 1, leaflet1[i], lipidNameList[i])
						orderParam.append(down)
					ordAvg = sum(orderParam) / len(orderParam)
					ordErr = statistics.pstdev(orderParam) / 2
					resultFile.write(f"Average order parameter (Lower) for {lipidNameList[i]} = {ordAvg} with a standard error of {ordErr}\n")	
				except:
					numfailures += 1

		# find lateral diffusion with gmx msd. Requires knowing head groups for lipids, which
		# is calculated below
		headGroups = []
		for i in range(len(lipidNameList)):
			lipid = lipidNameList[i]
			if leaflet1[i] > 0:
				if lipid[2:4] == "PC" or lipid[2:4] == "PE" or lipid[2:4] == "PS" or lipid[2:4] == "PG" or lipid[2:4] == "PA" or lipid[2:4] == "SM" or lipid[2:4] == "P6":
					headGroups.append(f"1 & r {lipid} & a PO4")
				elif lipid[2:4] == "PI":
					headGroups.append(f"1 & r {lipid} & a CP")
				elif lipid[2:4] == "P1":
					headGroups.append(f"1 & r {lipid} & a P1")
				elif lipid[2:4] == "P2":
					headGroups.append(f"1 & r {lipid} & a P2")
				elif lipid[2:4] == "P3":
					headGroups.append(f"1 & r {lipid} & a P3")
				elif lipid == "CHOL":
					headGroups.append(f"1 & r {lipid} & a ROH")
				else:
					print(f"error: {lipid} not of a type included in this script")
			if leaflet2[i] > 0:
				if lipid[2:4] == "PC" or lipid[2:4] == "PE" or lipid[2:4] == "PS" or lipid[2:4] == "PG" or lipid[2:4] == "PA" or lipid[2:4] == "SM" or lipid[2:4] == "P6":
					headGroups.append(f"0 & r {lipid} & a PO4")
				elif lipid[2:4] == "PI":
					headGroups.append(f"0 & r {lipid} & a CP")
				elif lipid[2:4] == "P1":
					headGroups.append(f"0 & r {lipid} & a P1")
				elif lipid[2:4] == "P2":
					headGroups.append(f"0 & r {lipid} & a P2")
				elif lipid[2:4] == "P3":
					headGroups.append(f"0 & r {lipid} & a P3")
				elif lipid == "CHOL":
					headGroups.append(f"0 & r {lipid} & a ROH")
				else:
					print(f"error: {lipid} not of a type included in this script")
		
		command = "(echo name 0 up; echo name 1 down; "		
		for group in headGroups:
			command += 'echo "' + group + '"; '
		command += "echo q) | " + f"gmx make_ndx -f {tprFP} -n {runDir}/analysis/bilayer_leaflet_0000.ndx -o {runDir}/analysis/diffusion.ndx"
		subprocess.call(command, shell=True)

		command = "("
		for i in range(len(headGroups)):
			command += f"echo {i+2}; "
		command = command[:-2] + f") | gmx msd -s {tprFP} -f {xtcFP} -lateral z -n {runDir}/analysis/diffusion.ndx -ngroup {len(headGroups)} -o {runDir}/analysis/diffusion_coeffs.xvg -b {finalSimulationTime / 5} -e {finalSimulationTime}"
		subprocess.call(command, shell=True)


		# find diffusion values from the output xvg files and write them to our output file
		diffusionFile = open(runDir + "/analysis/diffusion_coeffs.xvg", "r")
		for line in diffusionFile.readlines():
			if line[0:4] == "# D[":
				tokens = line[4:].split("_&_")
				layer = tokens[0]
				lipid = tokens[1]
				coeff = line.split("=")[1].split("(")[0]
				err = line.split("+/-")[1].split(")")[0]
				if layer == "up":
					resultFile.write("Lateral Diffusion Coef (Upper) for %s = %s with a standard error of %s\n" % (lipid, coeff, err))
				else:
					resultFile.write("Lateral Diffusion Coef (Lower) for %s = %s with a standard error of %s\n" % (lipid, coeff, err))

		# call cholesterol flip flop script
		try:
			numUp = 0
			numDown = 0
			subprocess.call(["python3", "trxcounter-bilayer-center_p2.py", "-f", xtcFP, "-s", tprFP, "-b", '0', "-e", str(finalSimulationTime), "-o", runDir + "/analysis/dist-chol-all.out"])
			cholFile = open(runDir + "/analysis/dist-chol-all.out", "r")
			line_counter = 0
			for line in cholFile.readlines():
				line_counter += 1
				if line_counter == 2:
					data = line.split("(")[1].split("down)")[0].split(" up/ ")
					numUp += int(data[0])
					numDown += int(data[1])
			cholFile.close()
			resultFile.write("CHOL %s flip-flops per ns (%s up/ %s down) \n\n" % (1000 * (numUp + numDown) / finalSimulationTime, numUp, numDown))
		except:
			numfailures += 1
		

		# call lateral pressure profile
		# Create a new .tpr file in the modified version of gromacs
		command = f"gmx_LS grompp -f {mdoFP} -c {groFP} -n {indFP} -p {topFP} -o {runDir}/20fs/20fs_LS.tpr"
		subprocess.call(command, shell=True)
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


                
	resultFile.close()

analyze([sys.argv[1]], "analysis/analysisResults" + sys.argv[1] + ".txt")

'''
# set up multiprocessing analysis

os.mkdir("analysis")
directoryList = os.listdir()
filteredList = list(filter((lambda x : x[0:3] == "run"), directoryList))

# run analysis in parallel with nprocs processors

nprocs = len(filteredList) 
ps = []
for i in range(nprocs):
	#a = int(i * len(filteredList) / nprocs)
	#b = int((i+1) * len(filteredList) / nprocs)
	p = multiprocessing.Process(target=analyze, args=([filteredList[i]], f"analysis/analysisResults{i+1}.txt"))
	ps.append(p)
	p.start()

for p in ps:
	p.join()

# clean up results

command = "cat "
for i in range(nprocs):
	command += f"analysis/analysisResults{i+1}.txt "
command += "> analysis/analysisResults.txt"
subprocess.call(command, shell=True)

for i in range(nprocs):
	subprocess.call(["rm", f"analysis/analysisResults{i+1}.txt"]) 
'''

