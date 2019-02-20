#!/usr/bin/python
# -*- coding: utf-8 -*-

# 2011.11.27 - Helgi I. Ingolfsson - Fix POPC and POPE (tail order is flipped in regular Martini 2.1 itp)
# 2015.07.01 - Alex H. de Vries - Adapt for lipidome (some lipid names changed, others not valid, 4 bead oleoyl tail)

from math import sqrt
from os import remove, system, path
from sys import argv, stdout
import subprocess

"""
  Compute (second rank) order parameter, defined as:

    P2 = 0.5*(3*<cos²(theta)> - 1)

  where "theta" is the angle between the bond and the bilayer normal.
  P2 = 1      perfect alignement with the bilayer normal
  P2 = -0.5   anti-alignement
  P2 = 0      random orientation

  Currently, the lipids DAPC, DLPC, DPPC, DOPC, POPC, DLPE, DPPE, DOPE, and POPE can be analyzed with this script.
  It is quite simple to extend the range by editing this program. 
  Search for DAPC and add similar lines for your lipid of interest. 
  Usage: %s <traj file> <tpr file> <initial time> <final time> <skip frames> <bilayer normal - xyz> <#lipids> <lipid type>

    > %s my.xtc my.tpr 5000 10000 5 0 1 0 64 DPPC

  will for example read the frames between 5 and 10 ns of the trajectory my.xtc, using information from my.tpr, expecting 64 DPPC lipids, calculating the order parameter for every 5th frame, averaging the results. 
  P2 will be calculated relative to the y-axis (the vector 0 1 0).

  The output is written to a file called order.dat.

  WARNING script will output all frames in one go, into files called frame_dump_XXX.gro and 
  then remove them so don't have any other files with this name in the current directory.
"""

### SOME TOOLS

# parse a .gro file
# return a list of coordinates
def read_gro(file, atoms):
  line_counter = 0
  number_of_particles = 0
  first, second = [], []
  for line in open(file):
    if line_counter == 1:
      number_of_particles = int(line)
    elif line_counter > 1 and line_counter < number_of_particles + 2:
      if line[10:15].strip() == atoms[0]:
        first.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
      elif line[10:15].strip() == atoms[1]:
        second.append([float(line[20:28]), float(line[28:36]), float(line[36:44])])
    line_counter += 1
  return [first, second]

### REAL STUFF

def computeOrderParameter(trajfile, tprfile, initial_time, final_time, traj_skip, normalx, normaly, normalz, number_of_lipids, lipid_type):


  # (normalized) orientation of bilayer normal
  orientation_of_bilayer_normal = [normalx, normaly, normalz]
  norm = sqrt(orientation_of_bilayer_normal[0]**2 + orientation_of_bilayer_normal[1]**2 + orientation_of_bilayer_normal[2]**2)
  for i in range(3):
    orientation_of_bilayer_normal[i] /= norm
  stdout.write("(Normalized) orientation of bilayer normal: ( %.3f | %.3f | %.3f ).\n" % (
    orientation_of_bilayer_normal[0], \
    orientation_of_bilayer_normal[1], \
    orientation_of_bilayer_normal[2]  \
  ))

  # output legend
  phosphatidylcholine_bond_names = " NC3-PO4 PO4-GL1 GL1-GL2 "
  phosphatidylethanolamine_bond_names = " NH3-PO4 PO4-GL1 GL1-GL2 "
  phosphatidylserine_bond_names = " CNO-PO4 PO4-GL1 GL1-GL2 "
  phosphatidylglycerol_bond_names = " GL0-PO4 PO4-GL1 GL1-GL2 "
  phosphatidic_acid_bond_names = " PO4-GL1 GL1-GL2 "
  phosphatidylinositol_bond_names = " C1-C2 C1-C3 C2-C3 C1-PO4 PO4-GL1 GL1-GL2 "
  # PCs
  if   lipid_type == "DAPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A GL2-D1B D1A-D2A D2A-D3A D3A-D4A D4A-C5A D1B-D2B D2B-D3B D3B-D4B D4B-C5B\n"
  elif lipid_type == "DBPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
  elif lipid_type == "DFPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A GL2-C1B C1B-D2B D2B-D3B D3B-D4B\n"
  elif lipid_type == "DGPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  elif lipid_type == "DIPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-D2B D2B-D3B D3B-C4B\n"
  elif lipid_type == "DLPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
  elif lipid_type == "DNPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"
  elif lipid_type == "DOPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-D2A D2A-C3A C3A-C4A C1B-D2B D2B-C3B C3B-C4B\n"
  elif lipid_type == "DPPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "DRPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-D6B\n"
  elif lipid_type == "DTPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A GL2-C1B C1B-C2B\n" 
  elif lipid_type == "DVPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-D3B D3B-C4B\n"
  elif lipid_type == "DXPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B\n"
  elif lipid_type == "DYPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-C3A GL2-C1B C1B-D2B D2B-C3B\n"
  elif lipid_type == "LPPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B\n"
  elif lipid_type == "PAPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PEPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PGPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PIPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "POPC": bond_names = phosphatidylcholine_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PRPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PUPC": bond_names = phosphatidylcholine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  # PEs
  elif lipid_type == "DAPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-C5B\n"
  elif lipid_type == "DBPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
  elif lipid_type == "DFPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A GL2-C1B C1B-D2B D2B-D3B D3B-D4B\n"
  elif lipid_type == "DGPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  elif lipid_type == "DIPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-D2B D2B-D3B D3B-C4B\n"
  elif lipid_type == "DLPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C1B-C2B C2B-C3B\n"
  elif lipid_type == "DNPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"
  elif lipid_type == "DOPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-D2A D2A-C3A C3A-C4A C1B-D2B D2B-C3B C3B-C4B\n"
  elif lipid_type == "DPPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A GL2-C1B C1A-C2A C2A-C3A C3A-C4A C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "DRPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-D6B\n"
  elif lipid_type == "DTPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A GL2-C1B C1B-C2B\n"
  elif lipid_type == "DUPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B\n"
  elif lipid_type == "DVPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-D3B D3B-C4B\n"
  elif lipid_type == "DXPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B\n"
  elif lipid_type == "DYPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-C3A GL2-C1B C1B-D2B D2B-C3B\n"
  elif lipid_type == "LPPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B\n"
  elif lipid_type == "PAPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PGPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PIPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "POPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PQPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PRPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PUPE": bond_names = phosphatidylethanolamine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  # PSs
  elif lipid_type == "DAPS": bond_names = phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-C5B\n"
  elif lipid_type == "DBPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
  elif lipid_type == "DFPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A GL2-C1B C1B-D2B D2B-D3B D3B-D4B\n"
  elif lipid_type == "DGPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  elif lipid_type == "DIPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-D2B D2B-D3B D3B-C4B\n"
  elif lipid_type == "DLPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A GL2-C1B C1B-C2B C2B-C3B\n"
  elif lipid_type == "DNPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"
  elif lipid_type == "DOPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-D2B D2B-C3B C3B-C4B\n"
  elif lipid_type == "DPPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "DRPS": bond_names = phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-D6B\n"
  elif lipid_type == "DTPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A GL2-C1B C1B-C2B\n"
  elif lipid_type == "DUPS": bond_names = phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B\n"
  elif lipid_type == "DVPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-D3B D3B-C4B\n"
  elif lipid_type == "DXPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B\n"
  elif lipid_type == "DYPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-C3A GL2-C1B C1B-D2B D2B-C3B\n"
  elif lipid_type == "LPPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B\n"
  elif lipid_type == "PAPS": bond_names = phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PGPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PIPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "POPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PQPS": bond_names = phosphatidylserine_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PRPS": bond_names = phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PUPS": bond_names = phosphatidylserine_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  # PGs
  elif lipid_type == "DAPG": bond_names = phosphatidylglycerol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-C5B\n"
  elif lipid_type == "DBPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
  elif lipid_type == "DFPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A GL2-C1B C1B-D2B D2B-D3B D3B-D4B\n"
  elif lipid_type == "DGPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  elif lipid_type == "DIPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-D2B D2B-D3B D3B-C4B\n"
  elif lipid_type == "DLPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A GL2-C1B C1B-C2B C2B-C3B\n"
  elif lipid_type == "DNPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"
  elif lipid_type == "DOPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-D2B D2B-C3B C3B-C4B\n"
  elif lipid_type == "DPPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "DRPG": bond_names = phosphatidylglycerol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-D6B\n"
  elif lipid_type == "DTPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A GL2-C1B C1B-C2B\n"
  elif lipid_type == "DVPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-D3B D3B-C4B\n"
  elif lipid_type == "DXPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B\n"
  elif lipid_type == "DYPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-C3A GL2-C1B C1B-D2B D2B-C3B\n"
  elif lipid_type == "LPPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B\n"
  elif lipid_type == "PAPG": bond_names = phosphatidylglycerol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PGPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PIPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "POPG": bond_names = phosphatidylglycerol_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PRPG": bond_names = phosphatidylglycerol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  # PAs
  elif lipid_type == "DAPA": bond_names = phosphatidic_acid_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-C5B\n"
  elif lipid_type == "DBPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B\n"
  elif lipid_type == "DFPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-D4A GL2-C1B C1B-D2B D2B-D3B D3B-D4B\n"
  elif lipid_type == "DGPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-D3B D3B-C4B C4B-C5B\n"
  elif lipid_type == "DIPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-D2B D2B-D3B D3B-C4B\n"
  elif lipid_type == "DLPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A GL2-C1B C1B-C2B C2B-C3B\n"
  elif lipid_type == "DNPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-D4A D4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-D4B D4B-C5B C5B-C6B\n"
  elif lipid_type == "DOPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-D2B D2B-C3B C3B-C4B\n"
  elif lipid_type == "DPPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "DRPA": bond_names = phosphatidic_acid_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-D1B D1B-D2B D2B-D3B D3B-D4B D4B-D5B D5B-D6B\n"
  elif lipid_type == "DTPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A GL2-C1B C1B-C2B\n"
  elif lipid_type == "DVPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-D3B D3B-C4B\n"
  elif lipid_type == "DXPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A C4A-C5A C5A-C6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B C4B-C5B C5B-C6B\n"
  elif lipid_type == "DYPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-C3A GL2-C1B C1B-D2B D2B-C3B\n"
  elif lipid_type == "LPPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B\n"
  elif lipid_type == "PAPA": bond_names = phosphatidic_acid_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PGPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-C2A C2A-D3A D3A-C4A C4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PIPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "POPA": bond_names = phosphatidic_acid_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PRPA": bond_names = phosphatidic_acid_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A D5A-D6A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PUPA": bond_names = phosphatidic_acid_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  # PIs
  elif lipid_type == "DPPI": bond_names = phosphatidylinositol_bond_names + "GL1-C1A C1A-C2A C2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PAPI": bond_names = phosphatidylinositol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PIPI": bond_names = phosphatidylinositol_bond_names + "GL1-C1A C1A-D2A D2A-D3A D3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "POPI": bond_names = phosphatidylinositol_bond_names + "GL1-C1A C1A-D2A D2A-C3A C3A-C4A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "PUPI": bond_names = phosphatidylinositol_bond_names + "GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-D5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
# Other Lipids
  elif lipid_type == "DPSM": bond_names = "NC3-PO4 PO4-AM1 AM1-AM2 AM1-T1A T1A-C2A C2A-C3A AM2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  elif lipid_type == "CHOL": bond_names = "C1-C2 ROH-R2 ROH-R3 R2-R3 C1-R2 C1-R3\n"
  elif lipid_type == "PAP6": bond_names = "C1-C2 C1-C3 C2-C3 C3-P4 C3-P5 C2-P5 C2-P4 C1-PO4 PO4-GL1 GL1-GL2 GL1-D1A D1A-D2A D2A-D3A D3A-D4A D4A-C5A GL2-C1B C1B-C2B C2B-C3B C3B-C4B\n"
  # output legend
  output_legend = "  Frame" + bond_names 

  # write the stuff
  stdout.write("\n " + output_legend)
  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n")
  output = open('order.dat', 'w')
  output.write(output_legend)
  output.write(("-"*(len(output_legend) - 1)) + "\n")

  # Output all frame using trjconv 
  stdout.write("Output all coordinate files \n")
  command = "echo %s | gmx trjconv -f %s -s %s -b %i -e %i -sep -skip %i -pbc whole -o frame_dump_.gro > /dev/null" % (lipid_type, trajfile, tprfile, initial_time, final_time, traj_skip)
  print(command)
  subprocess.call(command, shell=True)

  # For each dumped frame
  stdout.write("Starting P2 calculation")
  order_parameters = []
  file_count = 0
  bonds = []
  while True:
    filename = "frame_dump_" + str(file_count) + ".gro"
    if not path.isfile(filename) or path.getsize(filename) == 0:
        break
    
    stdout.write("Taking care of snapshot %s \n" % filename)

    # compute order parameter for each bond, for each snapshot
    current_order_parameters = []
    # bonds respectively involved in the head,
    #                             in the junction head-tail,
    #                             in each tail
    bonds = []

    for bond_name in bond_names.split():
      bonds.append(bond_name.split("-"))

    for bond in bonds:

      # parse .gro file, grep bead coordinates
      first, second = read_gro(filename, bond)

      # compute order parameter for each lipid
      order_parameter = 0.0
      for i in range(number_of_lipids):
        # vector between the two previous beads (orientation doesn't matter)
        vector = [0.0, 0.0, 0.0]
        for j in range(3):
          vector[j] = first[i][j] - second[i][j]
        norm2 = vector[0]**2 + vector[1]**2 + vector[2]**2
        # compute projection on the bilayer normal
        projection = vector[0]*orientation_of_bilayer_normal[0] + vector[1]*orientation_of_bilayer_normal[1] + vector[2]*orientation_of_bilayer_normal[2]
        # order parameter
        order_parameter += projection**2/norm2

      # compute final averaged order parameter
      # store everything in lists
      current_order_parameters.append(0.5*(3.0*(order_parameter/number_of_lipids) - 1.0))
    order_parameters.append(current_order_parameters)

    # write results
    results = "%7i" % file_count
    for order_parameter in current_order_parameters:
      results += "%8.3f" % order_parameter
    stdout.write(" " + results + "\n")
    output.write(results + "\n")

    remove(filename)
    file_count += 1
  # End while loop

  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")
  stdout.write("Snapshots analysis done.%s\n" % (" "*56))
  stdout.write("Computing averages...\n")

  # average order parameter
  averaged_order_parameters = []
  for i in range(len(bonds)):
    sum = 0.0
    for j in range(len(order_parameters)):
      sum += order_parameters[j][i]
    averaged_order_parameters.append(sum/len(order_parameters))
 
  # write results
  stdout.write("\n         " + bond_names)
  stdout.write(("-"*(len(output_legend) - 1)) + "\n")
  output.write(("-"*(len(output_legend) - 1)) + "\n")
  results = "average"
  for order_parameter in averaged_order_parameters:
    results += "%8.3f" % order_parameter
  stdout.write(" " + results + "\n")
  output.write(results + "\n")
  stdout.write(" " + ("-"*(len(output_legend) - 1)) + "\n\n")

  # Write abs average order parameters <Sn> (for carbon chains only)
  # WARNING this works with currenct lipids (all have defined x5 none carbon bonds) but for manually added lipids this might not be true
  ave_chain_s = 0
  for i in averaged_order_parameters[3:]: 
     ave_chain_s += abs(i)
  average_txt = "Abs average order parameters for carbon chains <Sn> = %8.3f \n\n" % (ave_chain_s / (len(averaged_order_parameters)-3))
  stdout.write(average_txt)
  output.write(average_txt)
  stdout.write("Results written in \"order.dat\".\n")

  output.close()

  return (ave_chain_s / (len(averaged_order_parameters)-3))
