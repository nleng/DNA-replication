#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys

from copy import copy
from math import sin,cos, pi
import numpy as np
from libreplihelpers import get_3d_points
import scipy.misc.pilutil as smp

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

if os.getenv('USER') == "garak":
	#vertfile = "results/paper_res_fork_0.txt"
	#beadfilename = "results/fullcalc14.txt"
	#beadfilename = "results/fullcalc15.txt"
	#beadfilename = "results/loops1.txt"
	vertfile = "/home/garak/phd/programs/replication/src/results/paper_hela_realistic_forks/res_fork_6000_280000_67000_28_500000_0_run.txt"
	beadfilename = "results/paper_random/res_res_1740.txt"
	
elif os.getenv('USER') == "loeb":
#	vertfile = "/home/loeb/phd/programs/replication/src/results/res_fork_0.0333333_400000_67000_5e+06_500000_0.txt"
	vertfile = "/home/loeb/phd/programs/replication/src/results/paper_hela_realistic/res_fork_6000_280000_67000_28_500000_0_run.txt"

	beadfilename = "/imports/hus.work/loeb/dna_metropolis/fullcalc15/res_res_9999.txt"
	#beadfilename = "/imports/hus.work/loeb/dna_metropolis/loops3/res_res_9999.txt"

#with_color = True
#onlydots = True
#radius = 5000.
#framenumber=350

paramslist = [	["D", 90, "early_S"],
		["E", 180, "mid_S"],
		["F", 351, "late_S"]
	]


for pent in [paramslist[-1]]:
	tmpdir = "results/images/" + pent[2]
	if not os.path.exists(tmpdir):
		os.makedirs(tmpdir)
	letter = pent[0]
	onlydots = True
	radius = 7500.
	
	framenumber= pent[1]
	
	minz = -3500
	#minz = 0
	maxz =  3500
	minx = -7500
	maxx =  7500
	deltaz = 150
	deltaxy = 40

	total_found = 0
	n_pixels = (maxx-minx)/deltaxy + 1

	print(float(maxx-minx)/deltaxy)
	print(float(maxz-minz)/deltaz)

	[gvertlist, rvertlist] = get_vertices_from_file(vertfile, framenumber)

	thelist = gvertlist + rvertlist
	
	print ("Started dryrun")
	# First, make a wasteful but necessary dry run to obtain the scaling factor.
	pvmax = 0
	for z in range(minz, maxz, deltaz):
		tmpdata = np.zeros( (n_pixels+40,n_pixels+40), dtype=np.uint16 )
		for entry in thelist:
			if entry[2] > z and entry[2] < z + deltaz:
				txpos = (entry[0]-minx)/deltaxy
				typos = (entry[1]-minx)/deltaxy
				tmpdata[20+txpos,20+typos] += 1
		
		scalefac = 0
		nnonzero = 0
		for i in range(0, n_pixels + 40):
			for j in range(0, n_pixels + 40):
				if tmpdata[i,j] > scalefac:
					scalefac = tmpdata[i,j]
				if not tmpdata[i,j] == 0:
					nnonzero +=1
		
		if scalefac == 0:
			scalefac = 1.
			
		if scalefac > pvmax:
			pvmax = scalefac
	
	print("Ended dryrun")
	listlist = [[gvertlist , 'green', 'eu'],[rvertlist, 'red', 'het']]

	for etr in listlist:
		cntr = 0
		for z in range(minz, maxz, deltaz):
			data = np.zeros( (n_pixels+40,n_pixels+40,3), dtype=np.uint8 )
			tmpdata = np.zeros( (n_pixels+40,n_pixels+40), dtype=np.uint16 )
			entrycounter = 0
			for entry in etr[0]:
				if entry[2] > z and entry[2] < z + deltaz:
					txpos = (entry[0]-minx)/deltaxy
					typos = (entry[1]-minx)/deltaxy
					tmpdata[20+txpos,20+typos] += 1
					entrycounter += 1
					total_found += 1
			print("Number of found entries", entrycounter)
			
			nnonzero = 0
			for i in range(0, n_pixels + 40):
				for j in range(0, n_pixels + 40):
					if not tmpdata[i,j] == 0:
						nnonzero +=1
			
			print("Nonzero pixels: ", nnonzero)
			
			for i in range(0, n_pixels + 40):
				for j in range(0, n_pixels + 40):
					tv = 255*(float(tmpdata[i,j])/float(pvmax))
					if etr[1] == 'green':
						data[i,j] = [0,tv,0]
					else:
						data[i,j] = [tv,0,0]
			
			ofname = "results/images/" + pent[2] + "/" + etr[2] + "_nr_" + str(int(cntr)) + ".tif"
			img = smp.toimage(data)
			#img.axes([minx,maxx,minx,maxx])	
			img.save(ofname)
			cntr += 1

	print(total_found)
