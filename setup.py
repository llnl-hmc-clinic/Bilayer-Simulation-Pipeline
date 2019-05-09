 #!/usr/bin/env python3
import os, configparser

dirName = ''

def main():
	config = configparser.ConfigParser()
	config.read('configuration.ini')
	dirName = os.path.expanduser(config['paths']['dir'])
	with open('files/header.txt', 'w') as f:
		f.write(f'#include "{dirName}/files/martini_v2.2.itp"\n')
		f.write(f'#include "{dirName}/files/martini_v2.0_ions.itp"\n')
		f.write(f'#include "{dirName}/files/martini_v2.0_all.itp"\n')
		f.write(f'#include "{dirName}/files/martini_v2.0_PAP6_02.itp"\n')
	os.chdir("..")
	parentDir = dirName = os.path.expanduser(config['paths']['parentDir'])
	folderName = dirName = os.path.expanduser(config['paths']['folderName'])
	os.system("cp -r 'Bilayer-Simulation-Pipeline' {0}".format(parentDir))

main()
