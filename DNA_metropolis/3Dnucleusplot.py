#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


from copy import copy
from math import sin,cos, pi
from libreplihelpers import get_3d_points

vertfile = "/home/corni/aaReplication/replicationSim/results/paper_hela_realistic_forks/res_fork_6000_280000_67000_28_500000_0_run.txt"
beadfilename = "/home/corni/aaReplication/dnaSim/outFile/res_res_1000.txt"	# gibts nicht!

#with_color = True
#onlydots = True
#radius = 5000.
#framenumber=350

paramslist = [	["D", 90, False, False],
		["E", 180, False, False],
		["F", 351, False, False],
		["G", 90, True, False],
		["H", 180, True, False],
		["I", 351, True, True]
	]

# hier zeit einstellen:
number = 5

background = 1.0
zellenrand = 0.0


letter = paramslist[number][0]
onlydots = True
#radius = 7500.
radius = 5000.	# hier 7500., die position der 2 mikrometer linie muss dann angepasst werden
#framenumber=351
framenumber= paramslist[number][1]
#framenumber=180
usecolor = paramslist[number][2]
plotlegend = paramslist[number][3]


def get_vertices_from_file(fname, toshow):
	infile = open(fname)

	counter = 0
	datalist = []
	for line in infile:
		if counter == toshow:
			zz = line.replace('\n','').split('\t')
			for i in range(1,len(zz)):
#				if float(zz[i]) < 225949719.0:
       				datalist.append(float(zz[i]))
			datalist.sort()
		counter += 1
	if not (os.path.exists(beadfilename)):
		print("MC results file not found.")
		sys.exit()

	tmplist = get_3d_points(datalist, beadfilename)

	print("done")
	return tmplist

def plotta(glista,rlista):
  eu = 0.
  hetero = 0.
  grenze = 10000	#1e308	# bei 1.79e+308 alle drin, ist das maximum von float
  xyz = [[], [], []]
  for i in range(len(glista)):
    for j in range(3):
      if glista[i][0]<grenze and glista[i][1]<grenze and glista[i][2]<grenze and glista[i][0]>-grenze and glista[i][1]>-grenze and glista[i][2]>-grenze:
	eu += 1./3.
	xyz[j].append(glista[i][j])

  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')
  ax.scatter(xyz[0], xyz[1], xyz[2], c='g', s=3, facecolor='g', lw = 0)
  
  xyz = [[], [], []]
  for i in range(len(rlista)):
    for j in range(3):
      if rlista[i][0]<grenze and rlista[i][1]<grenze and rlista[i][2]<grenze and rlista[i][0]>-grenze and rlista[i][1]>-grenze and rlista[i][2]>-grenze:
	hetero += 1./3.
	xyz[j].append(rlista[i][j])

  ax.scatter(xyz[0], xyz[1], xyz[2], c='r', s=3, facecolor='r', lw = 0)
  print eu
  print hetero
  print eu+hetero
  plt.show()
  
	
if __name__ == '__main__':
	global rvertlist
	global gvertlist

	[gvertlist, rvertlist] = get_vertices_from_file(vertfile, framenumber)
	plotta(gvertlist, rvertlist)
	print("N_eu:", len(gvertlist), "N_het:", len(rvertlist))

