"""
file:    generate.py
date:    2018 November 10
purpose: This program will serve the purpose of creating simulations of energy
		 minimization, 1fs, 5fs, and 15fs time step to relax the simulated lipid
		 bilayer for a 25fs or larger timestep simulations.

		 This program takes input from the file configurations.csv, which specifies
		 the lipid bilayers that will be simulated. This program will call 
		 GROMACS for producing the simulation results. configurations.csv has to
		 be edited as the given template, or else undefined behaviors might
		 occur.

         To launch a simulation, every line of the configurations file will be 
         read and used to generate the appropriate command line and arguments.

         Configuration format:
         using python3 to call generate.py will provide a template for configuration
         form.
"""

import sys,re,os,math, csv,random,subprocess,configparser,multiprocessing,tempfile,datetime, platform
import numpy as np

csvConfig = "./configurations.csv"#this is the input file in csv format
iniConfig = './configurations.ini'#this is the input file in ini format
data = {} #this collects the simulation parameters from the configuration file
next_run = 0#this keeps track of the number of runs we are doing
SYSTEM_SIZE = [25, 25, 15]#this is the system size of the simulation
run_type = sys.argv #this is the system arguments that either includes a run-type or not



""" importcsv imports a csv file that has specs for simulations and convert it 
	to a dictionary.
"""

def importcsv(csvfile):
	with open(csvfile, 'r') as fin:
		reader=csv.reader(fin, skipinitialspace=True, quotechar="'")
		for row in reader:
			if('simulation' in row[0]):
				data[row[0]]= {}
				for i in range(0,5):
					entry = next(reader)
					data[row[0]][entry[0]]= entry[1:]
					for j in range(0,len(entry)): #prevent empty cells from being read
						if(entry[j] == ''):
							data[row[0]][entry[0]]= entry[1:j]
							break

"""
	This provides a possibility of setting up simulations with ini configuration files.
"""

def importini(inifile):
	config = configparser.ConfigParser()
	config.read("inifile")
	for sec in config.sections(): #iterate through each simulation configuration
		for key, value in config.items(sec):
			data[sec] = {key: value.split(', ')}
	return data

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

"""
	Simulate is the function that produces and call commandlines from simulate.sh
	in order to run simulations.
"""

def simulate(lipids, bilayer):
	args = ""#arguments passed for initial simulation
	n = 0
	descr = ""#description of types of lipids that will be used in command lines
	counts = ""#counts for the number of lipids
	Utextnote = "upperLeaflet:"
	#a description text for upperleaflet in each folder that specifies the simulation
	Ltextnote = "lowerLeaflet:"
	#a description text for lowerleaflet in each folder that specifies the simulation
	Atextnote = "asymmetry: " + bilayer[2]
	i = 0
	for l in lipids: #create parameters for insane.py and gromacs.
		if int(bilayer[0][i]) > 0:
			args += "-u " + l + ":" + bilayer[0][i] + " "
			Utextnote += l + ": " + bilayer[0][i] + ", "
		if int(bilayer[1][i]) > 0:
			args += "-l " + l + ":" + bilayer[1][i] + " "
			Ltextnote += l + ": " + bilayer[1][i] + ", "
		if int(bilayer[0][i]) > 0 or bilayer[1][i] > 0:
			descr += l + " "
			counts += str(bilayer[0][i] + bilayer[1][i]) + " "
			n += 1
		i += 1
	args += "-pbc square -sol W -x {} -y {} -z {} -o bilayer.gro -p top.top ".format(SYSTEM_SIZE[0], SYSTEM_SIZE[1], SYSTEM_SIZE[2])
	args += "-asym " + str(bilayer[2])
	descr += "W "
	n += 1
	run_num = "run%04d" % run_number()
	if len(run_type) == 1: #provides options for doing relaxation simulations, 20fs simulations or both.
		if platform.system() == 'Darwin' or platform.system() == 'macosx':
			subprocess.call(["./generate-mac.sh", args, descr, str(n), str(bilayer[2]), counts, run_num, Utextnote, Ltextnote, Atextnote])
			subprocess.call(["./20fs-mac.sh", descr, counts, run_num, n, bilayer[3]])
		else:
			subprocess.call(["./generate.sh", args, descr, str(n), str(bilayer[2]), counts, run_num, Utextnote, Ltextnote, Atextnote])
			subprocess.call(["./20fs.sh", descr, counts, run_num, n, bilayer[3]])
	elif 'relax' in run_type:
		if platform.system() == 'Darwin' or platform.system() == 'macosx':
			subprocess.call(["./generate-mac.sh", args, descr, str(n), str(bilayer[2]), counts, run_num, Utextnote, Ltextnote, Atextnote])
		else:
			subprocess.call(["./generate.sh", args, descr, str(n), str(bilayer[2]), counts, run_num, Utextnote, Ltextnote, Atextnote])
	elif '20fs' in run_type:
		if platform.system() == 'Darwin' or platform.system() == 'macosx':
			subprocess.call(["./20fs-mac.sh", descr, counts, run_num, str(n), str(bilayer[3][0])])
		else:
			subprocess.call(["./20fs.sh", descr, counts, run_num, str(n), str(bilayer[3][0])])
"""
	In the main function we create a queue of simulations and run them through calling
	simulate().
"""

def main():
	importcsv(csvConfig)
	for key, value in data.items():
		lipidtype = value['lipid type']
		upper = value['upperLeaflet']
		lower = value['lowerLeaflet']
		asymmetry = value['asymmetry']
		steps = value['steps']
		total = 0.0
		for i in upper:
			total += float(i)
		queue = []
		if 'kc' in run_type and 'CHOL' in lipidtype: #check for when we want to keep the chlesterols
			CHOLindex = lipidtype.index('CHOL')
			CHOL = upper[CHOLindex]
			for asym in asymmetry:
				percentage = 1 - float(asym)/(total-int(CHOL))
				upper1 = []
				for i in range(len(upper)):
					if i == CHOLindex:
						upper1.append(str(upper[i]))
					else:
						upper1.append(str(round(float(upper[i])*(percentage))))
				queue.append([upper1, lower, asym, steps])
		else:
			for asym in asymmetry:
				percentage = 1 - float(asym)/total
				upper1 = []
				for i in range(len(upper)):
					upper1.append(str(round(float(upper[i])*(percentage))))
				queue.append([upper1, lower, asym, steps])
		while len(queue) > 0:
			e = queue.pop(0)
			simulate(lipidtype, e)

main()

