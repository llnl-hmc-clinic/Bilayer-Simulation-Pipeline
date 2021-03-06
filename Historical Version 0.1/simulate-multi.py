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

csvConfig = "./configurations.csv"#this is the input file in csv format
iniConfig = './configurations.ini'#this is the input file in ini format
data = {} #this collects the simulation parameters from the configuration file
initialRun = 0#this keeps track of the number of runs we are starting from
SYSTEM_SIZE = ["12", "12", "10"]#this is the system size of the simulation
runType = sys.argv #this is the system arguments that either includes a run-type or not
maxAttempt = 3 #this is the max number of attempts for simulations that failed
maxMulti = 4 #this is the max number of simulations that can be done simultaneously

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
	index of the current run and make changes to the global variable of nextRun.
"""
def runNumber():
	global nextRun
	n = nextRun
	nextRun += 1
	print(f"Set nextRun to {nextRun}")
	return n

"""
	Simulate is the function that call bash commands to manage the file 
	system and do gromacs simulations.
"""

def simulate(lipids, bilayer,runName):
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

	#call gromacs commands for energy minimization
	print("{0} starts".format(runName))
	os.system("mkdir {0}".format(runName))
	os.chdir(runName)
	tout = open('{0}.out'.format(runName), 'w')
	terr = open('{0}.err'.format(runName), 'w')
	tlog = open('{0}.log'.format(runName), 'w')
	tlog.write("{0} starts\n".format(runName))
	os.system("echo {0} >description.txt".format(Utextnote))
	os.system("echo {0} >> description.txt".format(Ltextnote))
	os.system("echo {0} >> description.txt".format(Atextnote))
	os.system("mkdir em")
	os.chdir("em")
	print("Start inserting membrane for {0}".format(runName))
	tlog.write("Start inserting membrane for {0}\n".format(runName))
	commands = "python2.7 ../../files/insane.py {0}".format(args)
	subprocess.call(commands, stdout=tout, stderr=terr, shell = True)
	if os.path.isfile("bilayer.gro"): 
		print("membrane inserted for {0}".format(runName))
		tlog.write("membrane inserted\n")
	else:
		print("membrane insertion failed for {0}".format(runName))
		tlog.write("membrane inserted\n")
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
	#call gromacs commands for the simulations following the energy minimization
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
	currentfile = "run{0}/{1}/{1}.xtc".format(n, dirname)
	if (not(os.path.isfile(currentfile))):
		os.system("rm -r {0}".format(runName))
		simulate(lipids, bilayer,runName)
	else:
		os.system("sbatch 20fs.sh 5 {0}/20fs".format(runName))

"""
	In the main function we create a queue of simulations and run them in parallel through 
	multithreading on simulate().
"""

def main():
	global initialRun	
	importcsv(csvConfig)
	queue = []
	#add the simulations into a queue
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
	#start the simulations from the queue
	ps = []
	while len(queue) > 0:
		e = queue.pop(0)
		for j in range(3):
			runName = "run" + str(n)
			p = multiprocessing.Process(target=simulate, args=(lipidtype, e, runName))
			ps.append(p)
			p.start()
			n += 1
			if len(ps) == maxMulti:
				for p in ps:
					p.join()
				ps = []
	for p in ps:
		p.join()

main()
