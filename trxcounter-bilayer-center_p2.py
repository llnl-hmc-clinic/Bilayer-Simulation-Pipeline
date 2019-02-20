#!/usr/bin/python
import numpy
import mdreader
import sys
import string
import math
from scipy.optimize import leastsq

# Manuels counter modified by Helgi 
# Helgi 2015.05.07 don't center on box middel but cacluate bilayer center using two gaussians - based on PO4 bead posisions

resultFile = open("dist-chol-all.out","w+")

# Simple gaussian function for curve fitting
def gaussian(x, mean, sd, a):
    norm = a * numpy.exp(-(x - mean)**2 / (2 * sd**2))
    return norm

# Fit x2 gaussians
def doublegfit(p, y, x):
    m1, m2, sd1, sd2, a1, a2 = p
    if m1 < 0 or m2 < 0 or m1 > 1000 or m2 > 1000 or a1 > 5 or a2 > 5 or a1 < 0 or a2 < 0:
        return 1e10
    else:
        y_fit = gaussian(x, m1, sd1, a1) + gaussian(x, m2, sd2, a2)
        return y - y_fit


syst = mdreader.MDreader()
syst.add_argument('-select', metavar='SELS', dest='selString', default='name ROH', help = 'string\tThe beads to select.')

# Bead has to be above or bellow this distance from center to be in on or the other bilayers
fliplim = 8

rohs = syst.select_atoms(syst.opts.selString)
allPO4 = syst.select_atoms("name PO4")

binRlenBinCounts = 800   
binRlenMax = 200 


# addition made by Michael 2-14-2019
i = int(len(allPO4.positions) / 2)
top = allPO4.positions[:i,2]
bot = allPO4.positions[i:,2]
topZ = float(sum(top) / len(top))
botZ = float(sum(bot) / len(bot))
topSD = float(math.sqrt(sum([pow(t - topZ, 2) for t in top]) / (len(top) - 1)))
botSD = float(math.sqrt(sum([pow(t - botZ, 2) for t in bot]) / (len(bot) - 1)))

plsq = numpy.array([botZ, topZ, botSD, topSD, 0.1, 0.1])     # Initial guesses for x2 gaussian used in leastsq (mean-inner, mean-outer, sd-inner, sd-outer, amplitude-inner, ampitude-outer) 

x = numpy.linspace(0,binRlenMax,binRlenBinCounts)

#for i in range(len(syst)):
def get_rohs_middle():
    global plsq
    cArray = allPO4.positions[:,2]
    y_raw,x_temp = numpy.histogram(cArray, bins=binRlenBinCounts, range=(0,binRlenMax)) # 0 is data, 1 is x spacing
    y = y_raw.astype(numpy.float32)/numpy.sum(y_raw)
    plsq = leastsq(doublegfit, plsq, args = (y,x))[0]
    return rohs.positions[:,2] - (plsq[1] + plsq[0]) / 2

cdx = numpy.array(syst.do_in_parallel(get_rohs_middle))

thresh = (cdx>fliplim).astype(numpy.int)
thresh -= cdx<-fliplim
numpy.savetxt("thresh.xvg", thresh, fmt="%d")
thresh = thresh.T

thresh = numpy.column_stack((thresh, numpy.ones(len(thresh))*1000))
thresh = thresh[thresh!=0] # Becomes 1D
diff = numpy.diff(thresh)

upflips = (diff==2).sum()
downflips = (diff==-2).sum()
numflip = upflips + downflips

totaltime = (syst.endframe - syst.startframe + 1) * syst.trajectory.dt / 1000.
res = "Counted %d moleculs/beads and %d frames over %.3f ns \n" % (len(rohs), syst.totalframes, totaltime)
resultFile.write(res)
res = "%.3f flip-flops per ns, from a total of %d events. (%d up/ %d down) \n" % (numflip/totaltime, numflip, upflips, downflips)
resultFile.write(res)

# Also count number in each leaflet and report 
thresh = (cdx>fliplim).astype(numpy.int)
countUpper = numpy.sum(thresh).astype(numpy.float)/syst.totalframes
thresh = (cdx<-fliplim).astype(numpy.int)
countLower = numpy.sum(thresh).astype(numpy.float)/syst.totalframes
res = "Of total of %d selected beads on average %.3f in upper and %.3f in lower leaflet \n" % (len(rohs), countUpper, countLower)
resultFile.write(res)

# Last line add summary 
res = "  -b  ,  -e  , n-beads , n-frames, numflips, countUpper, countLower, threshold \n"
resultFile.write(res)
res = "%.3f %.3f %d %d %d %.3f %.3f %.3f" % (syst.startframe*syst.trajectory.dt, syst.endframe*syst.trajectory.dt, len(rohs), syst.totalframes, numflip, countUpper, countLower, fliplim)
resultFile.write(res)

#res = "%.3f %d %d %d" % (numflip/totaltime, numflip, upflips, downflips)
#print(res)

resultFile.close()

# Print location of all beads with time 
with open(syst.opts.outfile, 'w') as OUT:
    OUT.write("#%s\n" %res)
    numpy.savetxt(OUT, numpy.column_stack((numpy.linspace(syst.startframe*syst.trajectory.dt, syst.endframe*syst.trajectory.dt, syst.totalframes), cdx/10)), fmt="%.3f")

