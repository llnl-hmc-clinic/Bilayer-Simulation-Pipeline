import sys,re,os,math,csv,random,subprocess,json,multiprocessing,tempfile,datetime
import numpy as np

CONFIGURATIONS = "./configurations.csv"#this is the output file

simulations = {'simulation1': 
				{'lipid type': ["POPC", "PAPC", "POPE", "DIPE", "DPSM", "PAPS", "PAP6", "CHOL"],
				'upperLeaflet': [243, 121, 20,  61, 242,   0,  0, 313],
				'lowerLeaflet': [139,  75, 54, 161, 108, 161, 22, 280],
				'asymmetry': [13, 23, 24],
				'sim0':['steep', '0.01', '10000', '2'],
				'sim1':['md', '0.001', '300000', '4'],
				'sim2':['md', '0.005', '300000', '6'],
				'sim3':['md', '0.015', '500000', '8'],
				'sim4':['md', '0.015', '500000', '8']},
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
