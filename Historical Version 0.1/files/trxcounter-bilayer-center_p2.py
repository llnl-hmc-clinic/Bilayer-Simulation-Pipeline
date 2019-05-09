#!/usr/bin/python
import numpy
import mdreader
import sys
import string
from scipy.optimize import leastsq

# Manuels counter modified by Helgi 
# Helgi 2015.05.07 don't center on box middel but cacluate bilayer center using two gaussians - based on PO4 bead posisions

# Simple gaussian function for curve fitting
def gaussian(x, mean, sd, a):
    norm = a * numpy.exp(-(x - mean)**2 / (2 * sd**2))
    return norm
    #norm = []
    #for i in range(len(x)):
    #    norm += [a * numpy.exp(-(x[i] - mean)**2 / (2 * sd**2))]
    #return numpy.array(norm)

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
#fliplim = 11
fliplim = 8

rohs = syst.select_atoms(syst.opts.selString)
allPO4 = syst.select_atoms("name PO4")
print "Running: " + string.join(map(str, sys.argv[:]))

#tseriesM = syst.timeseries(allPO4, props=["time","dimensions"],x=False, y=False)
#tseriesM.coords = tseriesM.coords.reshape(len(syst),len(allPO4)) 

#tseries = syst.timeseries(rohs, x=False, y=False)
#tseries.coords = tseries.coords.reshape(len(syst),len(rohs)) 

print("Fit bilayer middle")
binRlenBinCounts = 800   
binRlenMax = 200 
plsq = numpy.array([35., 70., 5., 5., 0.1, 0.1])     # Initial guesses for x2 gaussian used in leastsq (mean-inner, mean-outer, sd-inner, sd-outer, amplitude-inner, ampitude-outer) 
#plsq = numpy.array([40., 80., 5., 5., 0.1, 0.1])     # Initial guesses for x2 gaussian used in leastsq (mean-inner, mean-outer, sd-inner, sd-outer, amplitude-inner, ampitude-outer) 

#bilayerMiddle = numpy.zeros(len(syst))
x = numpy.linspace(0,binRlenMax,binRlenBinCounts)

#for i in range(len(syst)):
def get_rohs_middle():
    global plsq
    #cArray = tseriesM.coords[i,:]
    cArray = allPO4.positions[:,2]
    y_raw,x_temp = numpy.histogram(cArray, bins=binRlenBinCounts, range=(0,binRlenMax)) # 0 is data, 1 is x spacing
    y = y_raw.astype(numpy.float32)/numpy.sum(y_raw)
    plsq = leastsq(doublegfit, plsq, args = (y, x))[0]
    return rohs.positions[:,2] - (plsq[1] + plsq[0]) / 2

#bilayerMiddle = numpy.array(syst.do_in_parallel(get_bil_middle))
#cdx = tseries.coords - tseries.dimensions[:,2,None]/2
#cdx = tseries.coords - bilayerMiddle[:,None]
cdx = numpy.array(syst.do_in_parallel(get_rohs_middle))

thresh = (cdx>fliplim).astype(numpy.int)
thresh -= cdx<-fliplim
numpy.savetxt("thresh.xvg", thresh, fmt="%d")
thresh = thresh.T

#always_at_center = numpy.all(thresh==0, axis=0)
#if numpy.any(always_at_center):
#    raise ValueError("The following cholesterols never left the center:"
#                     "\n{}".format(numpy.where(always_at_center)[0]))
#while not numpy.all(thresh):
#    thresh[thresh==0] = numpy.roll(thresh, 1, axis=0)[thresh==0]

thresh = numpy.column_stack((thresh, numpy.ones(len(thresh))*1000))
thresh = thresh[thresh!=0] # Becomes 1D
diff = numpy.diff(thresh)

#thresh2 = numpy.roll(thresh, 1, axis=0)
#numflip = len(numpy.where((thresh-thresh2)[1:]!=0)[0])
upflips = (diff==2).sum()
downflips = (diff==-2).sum()
numflip = upflips + downflips

totaltime = (syst.endframe - syst.startframe + 1) * syst.trajectory.dt / 1000.
res = "Counted %d moleculs/beads and %d frames over %.3f ns" % (len(rohs), syst.totalframes, totaltime)
print res
res = "%.3f flip-flops per ns, from a total of %d events. (%d up/ %d down)" % (numflip/totaltime, numflip, upflips, downflips)
print res

# Also count number in each leaflet and report 
thresh = (cdx>fliplim).astype(numpy.int)
#countUpperArray = numpy.mean(thresh, axis=0) 
#countUpper = numpy.sum(countUpperArray)
countUpper = numpy.sum(thresh).astype(numpy.float)/syst.totalframes
thresh = (cdx<-fliplim).astype(numpy.int)
#countLowerArray = numpy.mean(thresh, axis=0)
#countLower = numpy.sum(countLowerArray)
countLower = numpy.sum(thresh).astype(numpy.float)/syst.totalframes
res = "Of total of %d selected beads on average %.3f in upper and %.3f in lower leaflet" % (len(rohs), countUpper, countLower)
print res

# Last line add summary 
res = "  -b  ,  -e  , n-beads , n-frames, numflips, countUpper, countLower, threshold "
print res
res = "%.3f %.3f %d %d %d %.3f %.3f %.3f" % (syst.startframe*syst.trajectory.dt, syst.endframe*syst.trajectory.dt, len(rohs), syst.totalframes, numflip, countUpper, countLower, fliplim)
print res

# Print location of all beads with time 
with open(syst.opts.outfile, 'w') as OUT:
    OUT.write("#%s\n" %res)
    numpy.savetxt(OUT, numpy.column_stack((numpy.linspace(syst.startframe*syst.trajectory.dt, syst.endframe*syst.trajectory.dt, syst.totalframes), cdx/10)), fmt="%.3f")

