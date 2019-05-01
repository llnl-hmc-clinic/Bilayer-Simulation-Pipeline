import sys,re,os,math,csv,random,subprocess,json,multiprocessing,tempfile,datetime
import numpy as np

CONFIGURATIONS = "./configurations.csv"#this is the output file
'''
The template for simulation is as the following:
In lipid type you specify the types of lipids that the simulated lipid  bilayer
has. Then in upperLeaflet and lowerLeaflet you specify the number of each type
of lipids in the same order. The asymmetry between two different leaflets for
each simulation is specified in asymmetry. From sim0 to simn you can specify the
type and length of simulations you use. The result of each simulation will be
used as the basis for the next simulation in order to relax the forces within the system. 
(i.e. The results of sim0 will be used as the basis for sim1.) The parameters are
integrator, time step, number of steps, and number of threads for the simulation.
The integrator is usually steep for energy minimization simulation, and md for the
following simulations. The time step has a unit of nanosecond. The number of thread
is up to the user, but suggested not to exceed the number of cores that the hardware
has.
'''
simulations = {'simulation1': 
				{'lipid type': ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"],
				'upperLeaflet': [243, 121, 20,  61, 242,   0,  0, 313],
				'lowerLeaflet': [139,  75, 54, 161, 108, 161, 22, 280],
				'asymmetry': [13, 23, 24],
				'sim0':['steep', '0.01', '10000', '2'],
				'sim1':['md', '0.001', '300000', '4'],
				'sim2':['md', '0.005', '300000', '6'],
				'sim3':['md', '0.015', '500000', '8'],
				'sim4':['md', '0.02', '500000', '8']},
		'simulation2': 
				{'lipid type': ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"],
				'lowerLeaflet': [243, 121, 20,  61, 242,   0,  0, 313],
				'upperLeaflet': [139,  75, 54, 161, 108, 161, 22, 280],
				'asymmetry': [13, 23, 24],
				'sim0':['steep', '0.01', '10000', '2'],
				'sim1':['md', '0.001', '300000', '4'],
				'sim2':['md', '0.005', '300000', '6'],
				'sim3':['md', '0.015', '500000', '8'],
                                'sim4':['md', '0.015', '500000', '8']},
		'end':{}
}



def main():
	with open(CONFIGURATIONS, 'w') as f:
		for key, value in simulations.items():
			f.write(str(key) + '\n')
			for key, value in value.items():
				f.write(str(key) + ',')
				for i in value:
					f.write(str(i) + ",")
				f.write('\n')
main()
