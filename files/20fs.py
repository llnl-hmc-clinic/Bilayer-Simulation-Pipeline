import sys,re,os,math,csv,random,subprocess,json,multiprocessing,tempfile,datetime
import numpy as np


CONFIGURATIONS = "./configurations.csv"#this is the input file
SIMULATION = "./20fs_and_analysis.sh"#this is the simulation file
data = {} #this collects the simulation parameters from the configuration file
next_run = 0#this keeps track of the number of runs we are doing

run_count = 0

def calculate_runs():
	existing = []
	for x in os.listdir('.'):
		m = re.search('run(\d+).*', x)
		if m:
			existing.append(int(m.group(1)))
	if existing:
		run_number = max(existing)

"""
    A helper function that keeps track of the number of runs. It returns the
    index of the current run.
"""
def run_number():
	global next_run
	n = next_run
	next_run += 1
	print(f"Set next_run to {next_run}")
	return n

""" importcsv imports a csv file that has specs for simulations and convert it 
	to a dictionary.
"""

def importcsv(csvfile):
	with open(csvfile, 'r') as fin:
		reader=csv.reader(fin, skipinitialspace=True, quotechar="'")
		for row in reader:
			space = 0
			for i in range(0,len(row)):
				if(row[i] == ''):
					data[row[0]]=row[1:i]
					space = 1
					break
				if(space == 0):
					data[row[0]]=row[1:]

"""
	Simulate is the function that produces and call commandlines from simulate.sh
	in order to run simulations.
"""

def simulate(lipids, bilayer,runnum):
	descr = ""#description of types of lipids that will be used in command lines
	counts = ""#counts for the number of lipids
	n = 0
	i = 0
	for l in lipids:
		if int(bilayer[0][i]) > 0 or bilayer[1][i] > 0:
			descr += l + " "
			counts += str(bilayer[0][i] + bilayer[1][i]) + " "
			n += 1
		i += 1
	descr += "W "
	run_num = "run%04d" % run_number()
	if platform.system() == 'Darwin' or 'macosx':
		subprocess.call(["./20fs-mac.sh", args, descr, str(n), str(bilayer[2]), counts, run_num, Utextnote, Ltextnote, Atextnote])
	else:
		subprocess.call(["./20fs.sh", args, descr, str(n), str(bilayer[2]), counts, run_num, Utextnote, Ltextnote, Atextnote])
"""
	In the main function we create a queue of simulations and run them through calling
	simulate().
"""

def main():
	global run_count
	importcsv(CONFIGURATIONS)
	numOfsim = len(data)//4
	for i in range(0, numOfsim):
		lipidtype = data['lipid type'+str(i+1)]
		upper = data['upperLeaflet'+str(i+1)]
		lower = data['lowerLeaflet'+str(i+1)]
		asymmetry = data['asymmetry'+str(i+1)]
		queue = []
		for asym in asymmetry:
			queue.append([upper, lower, asym])
		while len(queue) > 0:
			e = queue.pop(0)
			simulate(lipidtype, e, run_count)
			if not os.path.isdir(str(e)):
				queue.append(e)

main()
