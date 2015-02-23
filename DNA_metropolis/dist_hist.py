#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt

#infile = open("results.txt")
#squarelist = []
#for line in infile:
	#zz = line.replace('\n','').split('\t')
	#for k in range(0, len(zz)/3):
		#squarelist.append(float(zz[3*k])*float(zz[3*k]) + float(zz[3*k+1])*float(zz[3*k+1]))
	#break # only the first line

idir = "/imports/bott.work/loeb/dna_metropolis/nonucleus/"
for k in range(0,659):
	infile = open(idir + "res_res_" + str(k) + ".txt")
	xlist = []
	tlist = []
	for line in infile:
		zz = line.replace('\n','').split('\t')
		xpos = float(zz[0])
		ypos = float(zz[1])
		zpos = float(zz[2])
		for i in range(0, len(zz)/3-1):
			dx = float(zz[3*(i+1)])-xpos
			dy = float(zz[3*(i+1)+1])-ypos
			dz = float(zz[3*(i+1)+2])-zpos
			dist = sqrt(dx*dx+dy*dy+dz*dz)
			index = int(dist/10.)
			j = len(xlist)
			while index >= len(xlist):
				xlist.append(j*10)
				tlist.append(0)
				j += 1
			tlist[index] += 1
	
	infile.close()
	x = np.zeros(len(tlist))
	y = np.zeros(len(tlist))
	
	for i in range(0,len(tlist)):
		x[i] = xlist[i]
		y[i] = tlist[i]

	plt.clf()
	plt.plot(x,y,'rx')
	plt.axis([0,9000,0,300])
#plt.plot(x2,y2,'b')
#plt.xlabel('Genomic distance in Mbp')
#plt.ylabel('Mean square distance ($\mu m^2$)')
#plt.axis([0., 10., 0., 5.])
	plt.savefig(idir + "dist_hist_" + str(k) + ".png")
