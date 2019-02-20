import numpy as np 
import math
import os
import subprocess
import re
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})

#reading in all the property values into dictionaries
inputFile = open("analysis/analysisResults.txt", "r")
#orderdata = {}
areacompdata = {}
thicknessdata = {}
apldata = {}
aplerror = {}
flipflopdata = {}
for line in inputFile.readlines():
	if line[0:9] == "Asymmetry":
		asym = int(line[9:])
		'''
	if line[0:27] == "Average order parameter for":
		strList = line[27:].split("=")
		lipidName = strList[0][1:5]
		orderParamValue = float(strList[1])
		if lipidName not in orderdata:
			orderdata[lipidName] = {}
		if asym not in orderdata[lipidName]:
			orderdata[lipidName][asym] = []
		orderdata[lipidName][asym].append(orderParamValue)
		'''
	if line[0:23] == "Area compressibility = ":
		strList = line[23:].split("w")
		ACompValue = float(strList[0])
		if asym not in areacompdata:
			areacompdata[asym] = []
		areacompdata[asym].append(ACompValue)
	if line[0:20] == "Bilayer Thickness = ":
		strList = line[19:].split("w")
		ThicknessValue = float(strList[0])
		if asym not in thicknessdata:
			thicknessdata[asym] = []
		thicknessdata[asym].append(ThicknessValue)
	if line[0:17] == "Area per Lipid = ":
		strList = line[17:].split("w")
		AplValue = float(strList[0])
		if asym not in apldata:
			apldata[asym] = []
		apldata[asym].append(AplValue)
		AplError = float(strList[1][24:])
		if asym not in aplerror:
			aplerror[asym] = []
		aplerror[asym].append(AplError)
	if line[0:1].isdigit():
		strList = line.split("f")
		FlipFlopValue = float(strList[0])
		if asym not in flipflopdata:
			flipflopdata[asym] = []
		flipflopdata[asym].append(FlipFlopValue)


#getting the average order parameter for simulations of the same asymmetry
'''for lipid in orderdata:
	for asym in orderdata[lipid]:
		orderdata[lipid][asym] = sum(orderdata[lipid][asym])/len(orderdata[lipid][asym])
#make plot for order parameter
for lipid in orderdata:
	asymmetries = sorted(orderdata[lipid])
	orderparameters = []
	for asym in asymmetries:
		orderparameters.append(orderdata[lipid][asym])
	print(orderparameters)
	plt.plot(asymmetries, orderparameters, label=lipid, marker='o')'''

for asym in thicknessdata:
	thicknessdata[asym] = sum(thicknessdata[asym])/len(thicknessdata[asym])
thickasymmetries = sorted(thicknessdata)
thickness = []
for asym in thickasymmetries:
	thickness.append(thicknessdata[asym])

for asym in apldata:
	apldata[asym] = sum(apldata[asym])/len(apldata[asym])
aplasymmetries = sorted(apldata)
areaperlipid = []
for asym in aplasymmetries:
	areaperlipid.append(apldata[asym])

for asym in areacompdata:
	areacompdata[asym] = sum(areacompdata[asym])/len(areacompdata[asym])
acompasymmetries = sorted(areacompdata)
acomp = []
for asym in acompasymmetries:
	acomp.append(areacompdata[asym])

for asym in flipflopdata:
	flipflopdata[asym] = sum(flipflopdata[asym])/len(flipflopdata[asym])
ffasymmetries = sorted(flipflopdata)
flipflop = []
for asym in ffasymmetries:
	flipflop.append(flipflopdata[asym])

for asym in aplerror:
	toterr = 0
	for err in aplerror[asym]:
		toterr += err**2
	toterr /= len(aplerror[asym])
	aplerror[asym] = toterr
areaperlipiderror = []
for asym in aplasymmetries:
	areaperlipiderror.append(aplerror[asym])

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12,8))
ax1.plot(thickasymmetries, thickness, marker='o', markersize=4, color='tab:orange')
ax2.plot(aplasymmetries, areaperlipid, marker='o', markersize=4, color='tab:orange')
ax3.plot(acompasymmetries, acomp, marker='o', markersize=4, color='tab:orange')
ax4.plot(ffasymmetries, flipflop, marker='o', markersize=4, color='tab:orange')
ax1.set_xlabel("Asymmetry")
ax2.set_xlabel("Asymmetry")
ax3.set_xlabel("Asymmetry")
ax4.set_xlabel("Asymmetry")
#ax1.set_xticks([])
#ax2.set_xticks([])
#ax3.set_xticks([])
#ax4.set_xticks([])
#ax1.set_yticks([])
#ax2.set_yticks([])
#ax3.set_yticks([])
#ax4.set_yticks([])
ax1.set_ylabel("Bilayer Thickness")
ax2.set_ylabel("Area Per Lipid")
ax3.set_ylabel("Area Compressibility")
ax4.set_ylabel("Cholesterol Flip-Flop")
plt.savefig("multiplot.png")

plt.clf()

#plt.plot(aplasymmetries, areaperlipid, color='tab:orange', marker='o')
plt.errorbar(aplasymmetries, areaperlipid, yerr=areaperlipiderror, color='tab:orange')
plt.xlabel("Asymmetry")
plt.ylabel("Area per Lipid")
plt.savefig("detailplot.png")
print(areaperlipiderror)
