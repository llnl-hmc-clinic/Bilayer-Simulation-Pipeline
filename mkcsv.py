import sys,re,os,math,csv,random,subprocess,json,multiprocessing,tempfile,datetime
import numpy as np

CONFIGURATIONS = "./configurations.csv"#this is the output file

simulations = {'simulation1': 
				{'lipid type': ['DOPC', 'POPC', 'PAPE'],
				'upperLeaflet': ['300', '300', '400'],
				'lowerLeaflet': ['300', '300', '400'],
				'asymmetry': ['0', '25', '50', '75'],
				'steps': '15000000'},
				'simulation2': 
				{'lipid type': ['DPPC', 'POPC', 'PAPE'],
				'upperLeaflet': ['200', '300', '500'],
				'lowerLeaflet': ['200', '300', '500'],
				'asymmetry': ['0', '25', '50', '75', '100'],
				'steps': '15000000'}
}



def main():
	with open(CONFIGURATIONS, 'w') as f:
		for key, value in simulations.items():
			f.write(key + '\n')
			f.write('lipid type,')
			for lipid in value['lipid type']:
				f.write(lipid + ",")
			f.write('\nupperLeaflet,')
			for up in value['upperLeaflet']:
				f.write(str(up) + ",")
			f.write('\nlowerLeaflet,')
			for low in value['lowerLeaflet']:
				f.write(str(low) + ",")
			f.write('\nasymmetry,')
			for asym in value['asymmetry']:
				f.write(str(asym) + ",")
			f.write('\nsteps,')
			f.write(value['steps'] + ",")
			f.write('\n')
main()