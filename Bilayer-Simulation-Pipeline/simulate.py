"""
file:    generate.py
date:    2018 November 10
author:  HMC-LLNL-Clinic-Team
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

import sys,re,os,math, csv,random,subprocess,configparser,multiprocessing,tempfile,datetime, platform, analysis
import numpy as np

csvConfig = "./configurations.csv"#this is the input file in csv format
iniConfig = './configurations.ini'#this is the input file in ini format
data = {} #this collects the simulation parameters from the configuration file
next_run = 0#this keeps track of the number of runs we are doing
SYSTEM_SIZE = ["12", "12", "10"]#this is the system size of the simulation
run_type = sys.argv #this is the system arguments that either includes a run-type or not
maxAttempt = 3


""" isLast function is the function that checks whether this entry is the last 
	entry. If it is, it returns true as well as the last entry.
"""

def isLast(itr):
  old = itr.next()
  for new in itr:
    yield False, old
    old = new
  yield True, old

""" importcsv imports a csv file that has specs for simulations and convert it 
	to a dictionary.
"""

def importcsv(csvfile):
	with open(csvfile, 'r') as fin:
		reader=csv.reader(fin, skipinitialspace=True, quotechar="'")
		entry = [0]
		for row in reader:
			if('simulation' in row[0] or 'simulation' in entry[0]):				
				if('simulation' in row[0]):
					entry = next(reader)
				if('simulation' in entry[0]):
					exchange = row
					row = entry
					entry = exchange
				data[row[0]]= {}
				#keeps recording data in one simulation
				while 'simulation' not in entry[0] and 'end' not in entry[0]:
					data[row[0]][entry[0]]= entry[1:]
					for j in range(0,len(entry)): 
					#prevents empty cells from being read
						if(entry[j] == ''):
							data[row[0]][entry[0]]= entry[1:j]
							break
					entry = next(reader)

"""
	A helper function that keeps track of the number of runs. It returns the
	index of the current run and make changes to the global variable of next_run.
"""
def runNumber():
	global next_run
	n = next_run
	next_run += 1
	print(f"Set next_run to {next_run}")
	return n

"""
	Simulate is the function that call bash commands to manage the file 
	system and do gromacs simulations.
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
		if int(bilayer[0][i]) > 0 or int(bilayer[1][i]) > 0:
			descr += l + " "
			counts += str(bilayer[0][i] + str(bilayer[1][i])) + " "
			n += 1
		i += 1
	descr += "W "
	args += "-pbc square -sol W -salt 0.15 -x {} -y {} -z {} -o bilayer.gro -p top.top ".format(SYSTEM_SIZE[0], SYSTEM_SIZE[1], SYSTEM_SIZE[2])
	args += "-asym " + str(bilayer[2])
	n += 1
	run_num = "run" + str(runNumber())
	tout = open(f'{run_num}.out', 'w')
	terr = open(f'{run_num}.err', 'w')

	#call gromacs commands for simulations
	os.system("mkdir {0}".format(run_num))
	os.chdir(run_num)
	os.system("echo {0} >description.txt".format(Utextnote))
	os.system("echo {0} >> description.txt".format(Ltextnote))
	os.system("echo {0} >> description.txt".format(Atextnote))
	os.system("mkdir em")
	os.chdir("em")
	commands = "python2.7 ../../files/insane.py {0}".format(args)
	subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
	os.system("cat ../../files/header.txt top.top > topol.top")
	if platform.system() == 'Darwin' or platform.system() == 'macosx':	
		os.system("sed -i ' ' '5d' topol.top") #for mac
	else:		
		os.system("sed -i '5d' topol.top") #for linux 
	os.system("cp ../../files/martini_v2.x_new-rf.mdp ../em/")
	if platform.system() == 'Darwin' or platform.system() == 'macosx':
		os.system('sed -i " " "s/REPLACE1/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][0]))
		os.system('sed -i " " "s/REPLACE2/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][1]))
		os.system('sed -i " " "s/REPLACE3/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][2]))
	else:
		os.system('sed -i "s/REPLACE1/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][0]))
		os.system('sed -i "s/REPLACE2/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][1]))
		os.system('sed -i "s/REPLACE3/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][2]))
	commands = '(echo del 1-200; echo "r W | r NA+ | r CL-"; echo name 1 Solvent; echo !1; echo name 2 Membrane; echo q) | gmx make_ndx -f bilayer.gro -o index.ndx'
	subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
	commands = "gmx grompp -f martini_v2.x_new-rf.mdp -c bilayer.gro -p topol.top -n index.ndx -o em.tpr"
	subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
	commands = "gmx mdrun -deffnm em -v -nt {0} -dlb yes".format(bilayer[3][3])
	subprocess.call(commands, stdout=tout, stderr=terr, shell = True)	
	os.chdir("..")
	dirname = "em"
	for i in range(len(bilayer)-4):
		lastdir = dirname
		dirname = str(int(float(bilayer[i+4][1])*1000)) + "fs"
		os.system("mkdir {0}".format(dirname))#create a new directory named after the timestep of the simulation
		os.system("cp {0}/{0}.gro {0}/topol.top {1}".format(lastdir, dirname))	
		os.chdir(dirname)
		os.system("cp ../../files/martini_v2.x_new-rf.mdp ../{0}/index.ndx ../{1}/".format(lastdir, dirname))
		if platform.system() == 'Darwin' or platform.system() == 'macosx':	
			os.system('sed -i " " "s/REPLACE1/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][0]))
			os.system('sed -i " " "s/REPLACE2/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][1]))
			os.system('sed -i " " "s/REPLACE3/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][2]))
		else:
			os.system('sed -i "s/REPLACE1/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][0]))
			os.system('sed -i "s/REPLACE2/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][1]))
			os.system('sed -i "s/REPLACE3/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][2]))

		commands = "gmx grompp -f martini_v2.x_new-rf.mdp -c {0}.gro -p topol.top -n index.ndx -o {1}.tpr".format(lastdir, dirname)
		subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
		commands = "gmx mdrun -deffnm {0} -v -nt {1} -dlb yes".format(dirname, bilayer[i+4][3])
		subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
		os.chdir("..")
	os.chdir("..")

"""
	In the main function we create a queue of simulations and run them through calling
	simulate().
"""

def main():
	global next_run
	importcsv(csvConfig)
	existing = []
	for x in os.listdir('.'):
		m = re.search('run(\d+).*', x)
		if m:
			existing.append(int(m.group(1)))
	if existing:
		next_run = 1 + max(existing)
	else:
		next_run = 0
	queue = []
	for key, value in data.items():
		lipidtype = value['lipid type']
		upper = value['upperLeaflet']
		lower = value['lowerLeaflet']
		asymmetry = value['asymmetry']
		NoS = len(value)-4 #the number of keys
		for asym in asymmetry:
			simulation = [upper, lower, asym]
			for i in range(NoS):
				simulation.append(value['sim{0}'.format(i)])
			queue.append(simulation)
		while len(queue) > 0:
			e = queue.pop(0)
			for j in range(3):
				attempts = 0
				simulate(lipidtype, e)
				dirname = str(int(float(e[-1][1])*1000)) + "fs"
				currentfile = "run{0}/{1}/{1}.xtc".format(next_run-1, dirname)
				while (not(os.path.isfile(currentfile)) and attempts < 4):
					next_run -= 1
					simulate(lipidtype, e)
					attempts += 1
				if (attempts == 4):
					attempts = 0
				run_num = "run" + str(next_run - 1)
				os.system("python3 analysis.py {0}".format(run_num))

main()
