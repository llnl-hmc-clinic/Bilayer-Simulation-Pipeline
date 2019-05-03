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

import sys,re,os,math, csv,random,subprocess,configparser,multiprocessing,tempfile,datetime, platform
import numpy as np

filesDir = ''#this is the directory where we can get the needed files from
configFile = './configuration.ini'#this is the input file in form of configParser
data = {} #this collects the simulation parameters from the configuration file
initialRun = 0#this keeps track of the number of runs we are doing
systemSize = []#this is the system size of the simulation
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

"""importConfig() will import the configuration file and make changes to the
   global parameters that are configurized in the configuration file.
"""

def importConfig(iniFile):
	global systemSize
	global data
	global filesDir
	config = configparser.ConfigParser()
	config.read('configuration.ini')
	systemSize = config['general']['systemSize'].split(', ')
	filesDir = os.path.expanduser(config['paths']['dir'])
	for section in config:
		if (section[:10] == 'simulation'):
			data[section] = {}
			for name in config[section]:
				data[section][name] = config[section][name].split(', ')

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

def simulate(lipids, bilayer, runName):
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
	args += "-pbc square -sol W -salt 0.15 -x {} -y {} -z {} -o bilayer.gro -p top.top ".format(systemSize[0], systemSize[1], systemSize[2])
	args += "-asym " + str(bilayer[2])
	n += 1

	#call gromacs commands for simulations
	os.system("mkdir {0}".format(runName))
	os.chdir(runName)
	tout = open(f'{runName}.out', 'w')
	terr = open(f'{runName}.err', 'w')
	tlog = open('{0}.log'.format(runName), 'w')
	tlog.write("{0} starts\n".format(runName))
	os.system("echo {0} >description.txt".format(Utextnote))
	os.system("echo {0} >> description.txt".format(Ltextnote))
	os.system("echo {0} >> description.txt".format(Atextnote))
	os.system("mkdir em")
	os.chdir("em")
	print("Start inserting membrane for {0}".format(runName))
	tlog.write("Start inserting membrane for {0}\n".format(runName))
	commands = "python2.7 {1}/files/insane.py {0}".format(args, filesDir)
	subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
	if os.path.isfile("bilayer.gro"): 
		print("membrane inserted for {0}".format(runName))
		tlog.write("membrane inserted\n")
	else:
		print("membrane insertion failed for {0}".format(runName))
		tlog.write("membrane inserted\n")
	os.system("cat {0}/files/header.txt top.top > topol.top".format(filesDir))
	if platform.system() == 'Darwin' or platform.system() == 'macosx':	
		os.system("sed -i ' ' '5d' topol.top") #for mac
	else:
		os.system("sed -i '5d' topol.top") #for linux 
	os.system("cp {0}/files/martini_v2.x_new-rf.mdp ../em/".format(filesDir))
	if platform.system() == 'Darwin' or platform.system() == 'macosx':
		os.system('sed -i " " "s/REPLACE1/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][0]))
		os.system('sed -i " " "s/REPLACE2/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][1]))
		os.system('sed -i " " "s/REPLACE3/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][2]))
	else:
		os.system('sed -i "s/REPLACE1/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][0]))
		os.system('sed -i "s/REPLACE2/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][1]))
		os.system('sed -i "s/REPLACE3/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[3][2]))
	print("start energy minimization for {0}".format(runName))
	tlog.write("start energy minimization\n")
	commands = '(echo del 1-200; echo "r W | r NA+ | r CL-"; echo name 1 Solvent; echo !1; echo name 2 Membrane; echo q) | gmx make_ndx -f bilayer.gro -o index.ndx'
	subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
	commands = "gmx grompp -f martini_v2.x_new-rf.mdp -c bilayer.gro -p topol.top -n index.ndx -o em.tpr"
	subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
	commands = "gmx mdrun -deffnm em -v -nt {0} -dlb yes".format(bilayer[3][3])
	subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
	if os.path.isfile("em.gro"): 
		print("energy minimization finished for {0}".format(runName))
		tlog.write("energy minimization finished\n")
	else: 
		print("energy minimization failed for {0}".format(runName))
		tlog.write("energy minimization failed\n")	
	os.chdir("..")
	dirname = "em"
	for i in range(len(bilayer)-4):
		lastdir = dirname
		dirname = str(int(float(bilayer[i+4][1])*1000)) + "fs"
		os.system("mkdir {0}".format(dirname))#create a new directory named after the timestep of the simulation
		os.system("cp {0}/{0}.gro {0}/topol.top {1}".format(lastdir, dirname))	
		os.chdir(dirname)
		os.system("cp {2}/files/martini_v2.x_new-rf.mdp ../{0}/index.ndx ../{1}/".format(lastdir, dirname, filesDir))
		if platform.system() == 'Darwin' or platform.system() == 'macosx':	
			os.system('sed -i " " "s/REPLACE1/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][0]))
			os.system('sed -i " " "s/REPLACE2/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][1]))
			os.system('sed -i " " "s/REPLACE3/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][2]))
		else:
			os.system('sed -i "s/REPLACE1/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][0]))
			os.system('sed -i "s/REPLACE2/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][1]))
			os.system('sed -i "s/REPLACE3/{0}/" martini_v2.x_new-rf.mdp'.format(bilayer[i+4][2]))
		print("{0} simulation starts for {1}".format(dirname, runName))
		tlog.write("{0} simulation starts for {1}\n".format(dirname, runName))
		commands = "gmx grompp -f martini_v2.x_new-rf.mdp -c {0}.gro -p topol.top -n index.ndx -o {1}.tpr".format(lastdir, dirname)
		subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
		commands = "gmx mdrun -deffnm {0} -v -nt {1} -dlb yes".format(dirname, bilayer[i+4][3])
		subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
		if(os.path.isfile("{0}.gro".format(dirname))): 
			print("{0} simulation finished for {1}".format(dirname, runName))
			tlog.write("{0} simulation finished for {1}\n".format(dirname, runName))
		else: 
			print("{0} simulation failed for {1}".format(dirname, runName))
			tlog.write("{0} simulation failed for {1}\n".format(dirname, runName))
		os.chdir("..")
	os.chdir("..")
	dirname = str(int(float(bilayer[-1][1])*1000)) + "fs"
	currentfile = "{0}/{1}/{1}.xtc".format(runName, dirname)
	if (not(os.path.isfile(currentfile))):
		os.system("rm -r {0}".format(runName))
		simulate(lipids, bilayer,runName)
	else:
		currLocation = os.getcwd() 
		targetRun = currLocation + "/" + runName
		os.system("python3 {2}/analysis.py {0} {1}".format(runName, targetRun, filesDir))

"""
	In the main function we create a queue of simulations and run them through calling
	simulate().
"""

def main():
	global next_run
	importConfig(configFile)
	#check for the current folders and start new runs following those folders
	existing = []
	for x in os.listdir('.'):
		m = re.search('run(\d+).*', x)
		if m:
			existing.append(int(m.group(1)))
	if existing:
		n = max(existing) + 1
	else:
		n = 0
	queue = []
	#add the simulations into a queue
	for key, value in data.items():
		lipidtype = value['lipids']
		upper = value['upper']
		lower = value['lower']
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
				runName = "run" + str(n)
				simulate(lipidtype, e, runName)
	
main()
