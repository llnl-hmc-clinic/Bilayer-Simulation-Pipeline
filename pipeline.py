import sys,os,math,random,subprocess
import numpy as np

lipids=["DPPC", "DOPC", "POPC"]
exp1 = [[[1, 0, 0], [1, 0, 0], [0, 10, 20, 30, 40, 50]],
		[[0, 1, 0], [0, 1, 0], [0, 10, 20, 30, 40, 50]],
		[[0, 0, 1], [0, 0, 1], [0, 10, 20, 30, 40, 50]]]

def simulate(bilayer):
	args = ""
	n = 0
	descr = ""
	i = 0
	for l in lipids:
		if bilayer[0][i] > 0:
			args += "-u " + l + ":" + str(bilayer[0][i]) + " "
		if bilayer[1][i] > 0:
			args += "-l " + l + ":" + str(bilayer[1][i]) + " "
		if bilayer[0][i] > 0 or bilayer[1][i] > 0:
			descr += l + " "
			n += 1
		i += 1
	args += "-pbc square -sol W -x 7.5 -y 7.5 -z 10 -o bilayer.gro -p top.top "
	args += "-asym " + str(bilayer[2])
	descr += "W"
	n += 1
	subprocess.call(["./simulate.sh", args, descr, str(n), str(bilayer[2])])
	return

for e in exp1:
	for asym in e[2]:
		simulate([e[0], e[1], asym])

sys.exit(0)
