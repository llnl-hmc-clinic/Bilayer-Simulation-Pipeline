import configparser

config = configparser.ConfigParser()
config['simulation1'] = {'lipids': 'DOPC, POPC, PAPE',
                     'upperleaf': '300, 300, 400',
                     'lowerleaf': '300, 300, 400',
                     'asymmetry': '0, 25, 50, 75',
                     'steps': '1500000'}
config['simulation2'] = {'lipids': 'DOPC, POPC, PAPE',
                     'upperleaf': '200, 400, 400',
                     'lowerleaf': '200, 400, 400',
                     'asymmetry': '0, 25, 50, 75',
                     'steps': '1500000'}
with open('configuration.ini', 'w') as configfile:
	config.write(configfile)