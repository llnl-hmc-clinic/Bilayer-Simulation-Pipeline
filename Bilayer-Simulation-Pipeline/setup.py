import os, configparser

dirName = ''

def main():
	config = configparser.ConfigParser()
	config.read('configuration.ini')
	dirName = config['paths']['dir']
	os.chdir("..")
	retVal = os.getcwd()
	parentDir = config['paths']['parentDir']
	folderName = config['paths']['folderName']
	os.chdir(os.path.expanduser(parentDir))
	os.system("mkdir {0}".format(folderName))
	os.chdir(retVal)
	os.system("cp -r 'Bilayer-Simulation-Pipeline' {0}".format(dirName))

main()