import numpy as np 
import math
import os
import subprocess
import re
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

# reading in all the property values into dictionaries

inputFile = open("analysis/analysisResults.txt", "r")

btDict = {} # bilayer thickness
aplDict = {} # area per lipid
acompDict = {} # area compressibility
ordDict = {} # order parameters
diffDict = {} # lateral diffusion coefficient
ffDict = {} # cholesterol flip flop
'''
for lipid in lipidNameList:
	ordDict[lipid] = {}
	diffDict[lipid] = {}


for asym in asymList:
	btDict[asym] = []
	aplDict[asym] = {}
	aplDict[asym][0] = []
	aplDict[asym][1] = []
	acompDict[asym] = {}
	acompDict[asym][0] = []
	acompDict[asym][1] = []
	for lipid in lipidNameList:
		ordDict[lipid][asym] = {}
		ordDict[lipid][asym][0] = []
		ordDict[lipid][asym][1] = []
		diffDict[lipid][asym] = {}
		diffDict[lipid][asym][0] = []
		diffDict[lipid][asym][1] = []
	ffDict[asym] = []
'''

for line in inputFile.readlines():

	# on lines that list an asymmetry, switch asym to the new value
	if line.startswith("Asymmetry"):
		asym = int(line[9:])

	# parse bilayer thickness lines
	if line.startswith("Bilayer Thickness = "):
		strList = line[20:].split()
		value = float(strList[0])
		error = float(strList[6])
		if asym not in btDict:
			btDict[asym] = []
		btDict[asym].append((value, error))

	# parse area per lipid lines - by leaflet
	if line.startswith("Area per Lipid (Upper) = "):
		strList = line[25:].split()
		value = float(strList[0])
		error = float(strList[6])
		if asym not in aplDict:
			aplDict[asym] = {}
			aplDict[asym][0] = []
			aplDict[asym][1] = []
		aplDict[asym][0].append((value, error)) 	
	if line.startswith("Area per Lipid (Lower) = "):
		strList = line[25:].split()
		value = float(strList[0])
		error = float(strList[6])
		if asym not in aplDict:
			aplDict[asym] = {}
			aplDict[asym][0] = []
			aplDict[asym][1] = []
		aplDict[asym][1].append((value, error)) 

	# parse area compressibility lines - by leaflet
	if line.startswith("Area compressibility (Upper) = "):
		strList = line[31:].split()
		value = float(strList[0])
		error = float(strList[6])
		if asym not in acompDict:
			acompDict[asym] = {}
			acompDict[asym][0] = []
			acompDict[asym][1] = []
		acompDict[asym][0].append((value, error)) 
	if line[:31] == "Area compressibility (Lower) = ":
		strList = line[31:].split()
		value = float(strList[0])
		error = float(strList[6])
		if asym not in acompDict:
			acompDict[asym] = {}
			acompDict[asym][0] = []
			acompDict[asym][1] = []
		acompDict[asym][1].append((value, error))

	# parse order parameter lines - by leaflet and lipid
	if line.startswith("Average order parameter (Upper) for "):
		strList = line[36:].split()
		lipid = strList[0]
		value = float(strList[2])
		error = float(strList[8])
		if lipid not in ordDict:
			ordDict[lipid] = {}
		if asym not in ordDict[lipid]:
			ordDict[lipid][asym] = {}
			ordDict[lipid][asym][0] = []
			ordDict[lipid][asym][1] = []
		ordDict[lipid][asym][0].append((value, error))
	if line.startswith("Average order parameter (Lower) for "):
		strList = line[36:].split()
		lipid = strList[0]
		value = float(strList[2])
		error = float(strList[8])
		if lipid not in ordDict:
			ordDict[lipid] = {}
		if asym not in ordDict[lipid]:
			ordDict[lipid][asym] = {}
			ordDict[lipid][asym][0] = []
			ordDict[lipid][asym][1] = []
		ordDict[lipid][asym][1].append((value, error))

	# parse lateral diffusion lines - by leaflet and lipid
	if line.startswith("Lateral Diffusion Coef (Upper) for "):
		strList = line[35:].split()
		lipid = strList[0]
		value = float(strList[2])
		error = float(strList[8])
		if lipid not in diffDict:
			diffDict[lipid] = {}
		if asym not in diffDict[lipid]:
			diffDict[lipid][asym] = {}
			diffDict[lipid][asym][0] = []
			diffDict[lipid][asym][1] = []
		diffDict[lipid][asym][0].append((value, error))
	if line.startswith("Lateral Diffusion Coef (Lower) for "):
		strList = line[35:].split()
		lipid = strList[0]
		value = float(strList[2])
		error = float(strList[8])
		if lipid not in diffDict:
			diffDict[lipid] = {}
		if asym not in diffDict[lipid]:
			diffDict[lipid][asym] = {}
			diffDict[lipid][asym][0] = []
			diffDict[lipid][asym][1] = []
		diffDict[lipid][asym][1].append((value, error))

	# parse cholesterol flip flop lines
	if line.startswith("CHOL"):
		rate = float(line.split()[1])
		strList = line.split("flip-flops per ns (")[1].split()
		up = float(strList[0])
		down = float(strList[2])
		if asym not in ffDict:
			ffDict[asym] = []
		ffDict[asym].append((rate, up - down))

#plot thickness

btx = []
thickness = []
bterr = []
for asym in sorted(btDict):
	btx.append(asym)
	avgs = [t[0] for t in btDict[asym]]
	value = sum(avgs) / len(avgs)
	thickness.append(value)

	errs = [t[1] for t in btDict[asym]]
	e = sum(errs) / len(errs)
	bterr.append(e)

plt.figure()
plt.errorbar(btx, thickness, xerr=0, yerr=bterr, fmt='o')
plt.title("Bilayer Thickness (nm)")
#plt.show()
plt.savefig('graphs/thickness.png')

#plot apl

aplx = []
topapl = []
botapl = []
topaplerr = []
botaplerr = []
for asym in sorted(aplDict):
	aplx.append(asym)

	topavgs = [t[0] for t in aplDict[asym][0]]
	topvalue = sum(topavgs) / len(topavgs)
	topapl.append(topvalue)

	toperrs = [t[1] for t in aplDict[asym][0]]
	e = sum(toperrs) / len(toperrs)
	topaplerr.append(e)

	botavgs = [t[0] for t in aplDict[asym][1]]
	botvalue = sum(botavgs) / len(botavgs)
	botapl.append(botvalue)

	boterrs = [t[1] for t in aplDict[asym][1]]
	e = sum(boterrs) / len(boterrs)
	botaplerr.append(e)

plt.figure()
plt.errorbar(aplx, botapl, xerr=0, yerr=botaplerr, fmt='ro')
plt.errorbar(aplx, topapl, xerr=0, yerr=topaplerr, fmt='bo')
plt.title("Area Per Lipid")
#plt.show()
plt.savefig('graphs/apl.png')

#plot area compressibility

acompx = []
topareacomp = []
botareacomp = []
topacomperr = []
botacomperr = []
for asym in sorted(acompDict):
	acompx.append(asym)

	topavgs = [t[0] for t in acompDict[asym][0]]
	topvalue = sum(topavgs) / len(topavgs)
	topareacomp.append(topvalue)

	toperrs = [t[1] for t in acompDict[asym][0]]
	e = sum(toperrs) / len(toperrs)
	topacomperr.append(e)

	botavgs = [t[0] for t in acompDict[asym][1]]
	botvalue = sum(botavgs) / len(botavgs)
	botareacomp.append(botvalue)

	boterrs = [t[1] for t in acompDict[asym][1]]
	e = sum(boterrs) / len(boterrs)
	botacomperr.append(e)

plt.figure()
plt.errorbar(acompx, botareacomp, xerr=0, yerr=botacomperr, fmt='ro')
plt.errorbar(acompx, topareacomp, xerr=0, yerr=topacomperr, fmt='bo')
plt.title("Area Compressibility")
#plt.show()
plt.savefig('graphs/area_comp.png')

#plot order parameters

ordx = []
toporder = {}
botorder = {}
toporderr = {}
botorderr = {}
for lipid in ordDict:
	toporder[lipid] = []
	botorder[lipid] = []
	toporderr[lipid] = []
	botorderr[lipid] = []

print(ordDict.keys())
if len(ordDict.keys()) != 0:
	for asym in sorted(ordDict[list(ordDict.keys())[0]]):
		ordx.append(asym)

for lipid in ordDict:
	for asym in sorted(ordDict[lipid]):
		if len(ordDict[lipid][asym][0]) > 0:

			topavgs = [t[0] for t in ordDict[lipid][asym][0]]
			topvalue = sum(topavgs) / len(topavgs)
			toporder[lipid].append(topvalue)

			toperrs = [t[1]**2 for t in ordDict[lipid][asym][0]]
			e = math.sqrt(sum(toperrs)) / len(toperrs)
			toporderr[lipid].append(e)

		if len(orderdata[lipid][asym][1]) > 0:

			botavgs = [t[0] for t in ordDict[lipid][asym][1]]
			botvalue = sum(botavgs) / len(botavgs)
			botorder[lipid].append(botvalue)

			boterrs = [t[1]**2 for t in ordDict[lipid][asym][1]]
			e = math.sqrt(sum(boterrs)) / len(boterrs)
			botorderr[lipid].append(e)
# print(toporder["CHOL"])
# print(botorder["CHOL"])

for lipid in ordDict:
	if len(botorder[lipid]) > 0:
		plt.figure()
		plt.errorbar(ordx, botorder[lipid], xerr=0, yerr=botorderr[lipid], fmt='ro')
		plt.title("Order Parameter (Lower) for " + lipid)
		#plt.show()
		plt.savefig('graphs/order_bot_' + lipid + '.png')
	if len(toporder[lipid]) > 0 and not lipid == "CHOL":
		plt.figure()
		plt.errorbar(ordx, toporder[lipid], xerr=0, yerr=toporderr[lipid], fmt='bo')
		plt.title("Order Parameter (Upper) for " + lipid)
		#plt.show()
		plt.savefig('graphs/order_top_' + lipid + '.png')

#plot diffusion coefficients

diffx = []
topdiffusion = {}
botdiffusion = {}
topdifferr = {}
botdifferr = {}
for lipid in diffDict:
	topdiffusion[lipid] = []
	botdiffusion[lipid] = []
	topdifferr[lipid] = []
	botdifferr[lipid] = []

if len(diffDict.keys()) != 0:
	for asym in sorted(diffDict[list(diffDict.keys())[0]]):
		diffx.append(asym)

for lipid in diffDict:
	for asym in sorted(diffDict[lipid]):
		if len(diffDict[lipid][asym][0]) > 0:

			topavgs = [t[0] for t in diffDict[lipid][asym][0]]
			topvalue = sum(topavgs) / len(topavgs)
			topdiffusion[lipid].append(topvalue)

			toperrs = [t[1]**2 for t in diffDict[lipid][asym][0]]
			e = math.sqrt(sum(toperrs)) / len(toperrs)
			topdifferr[lipid].append(e)
		else:
			topdiffusion[lipid].append(0)
			topdifferr[lipid].append(0)


		if len(diffDict[lipid][asym][1]) > 0:

			botavgs = [t[0] for t in diffDict[lipid][asym][1]]
			botvalue = sum(botavgs) / len(botavgs)
			botdiffusion[lipid].append(botvalue)

			boterrs = [t[1]**2 for t in diffDict[lipid][asym][1]]
			e = math.sqrt(sum(boterrs)) / len(boterrs)
			botdifferr[lipid].append(e)
		else:
			botdiffusion[lipid].append(0)
			botdifferr[lipid].append(0)

for lipid in diffDict:
	if len(botdiffusion[lipid]) > 0:
		plt.figure()
		plt.errorbar(diffx, botdiffusion[lipid], xerr=0, yerr=botdifferr[lipid], fmt='ro')
		plt.title("Diffusion Coeffiecent (Lower) for " + lipid)
		#plt.show()
		plt.savefig('graphs/diffusion_bot_' + lipid + '.png')
	if len(topdiffusion[lipid]) > 0:
		plt.figure()
		plt.errorbar(diffx, topdiffusion[lipid], xerr=0, yerr=topdifferr[lipid], fmt='bo')
		plt.title("Diffusion Coeffiecent (Lower) for " + lipid)
		plt.savefig('graphs/diffusion_top_' + lipid + '.png')
		#plt.show()

#plot flip flop rate

ffx = []
flipflop = []
bias = []
for asym in sorted(ffDict):
	ffx.append(asym)
	avgs = [t[0] for t in ffDict[asym]]
	value = sum(avgs) / len(avgs)
	flipflop.append(value)

	bs = [t[1] for t in ffDict[asym]]
	b = sum(bs) / len(bs)
	bias.append(b)

plt.figure()
plt.plot(ffx, flipflop, 'o')
#plt.plot(x, bias, 'o')
plt.title("Flip flop rate")
#plt.show()
plt.savefig('graphs/flip.png')


# create multiplot figure, fill in graphs from above
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8), (ax9, ax10)) = plt.subplots(5, 2, figsize=(20,50))
# list of subplots to fill in
axList = (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10)
# parallel list of titles
propList = ("Order Parameters (Upper Leaflet)", "Order Parameter (Lower Leaflet)",
	"Lateral Diffusion Coefficient (Upper Leaflet)", "Lateral Diffusion Coefficient (Lower Leaflet)",
	"Area Compressibility (Upper Leaflet)", "Area Compressibility (Lower Leaflet)",
	"Area per Lipid (Upper Leaflet)", "Area per Lipid (Lower Leaflet)",
	"Bilayer Thickness", "Cholesterol Flip Flop")
# parallel lists of asym, data and error locations
asymList = (ordx, ordx, diffx, diffx, acompx, acompx, aplx, aplx, btx, ffx)
dataList = (toporder, botorder, topdiffusion, botdiffusion, topareacomp, botareacomp,
	topapl, botapl, thickness, flipflop)
errList = (toporderr, botorderr, topdifferr, botdifferr, topacomperr, botacomperr, 
	topaplerr, botaplerr, bterr, [])
# iterate through plots and fill in data
for ax, prop, x, data, err in zip(axList, propList, asymList, dataList, errList):
	ax.set_title(prop)
	if type(data) == dict:
		for lipid in data:
			ax.errorbar(x, data[lipid], xerr=0, yerr=err[lipid], fmt='ro')
	elif err == []:
		ax.plot(x, data, 'o')
	else:
		ax.errorbar(x, data, xerr=0, yerr=err, fmt='ro')
plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9, wspace=0.2, hspace=0.2)
plt.savefig("graphs/all_properties.png")