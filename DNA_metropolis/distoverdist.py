#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#infile = open("results.txt")
#squarelist = []
#for line in infile:
	#zz = line.replace('\n','').split('\t')
	#for k in range(0, len(zz)/3):
		#squarelist.append(float(zz[3*k])*float(zz[3*k]) + float(zz[3*k+1])*float(zz[3*k+1]))
	#break # only the first line

#infile = open("/imports/bott.work/loeb/dna_metropolis/startbig/res_distances_31.txt")
# infile = open("/imports/bott.work/loeb/dna_metropolis/res_distances_618.txt")
#infile = open("/imports/bott.work/loeb/dna_metropolis/newkappa/res_distances_543.txt")
infile = open("/imports/lydia.work/loeb/dna_metropolis/benchmark6_c4/res_distances_2520.txt")

tlist = []
for line in infile:
	zz = line.replace('\n','').split('\t')
	tlist.append([float(zz[0]), float(zz[1])])

tlist.sort()

x = np.zeros(len(tlist))
y = np.zeros(len(tlist))
x2 = np.zeros(len(tlist))
y2 = np.zeros(len(tlist))

for i in range(0,len(tlist)):
	x[i] = tlist[i][0]/1.e6
	y[i] = tlist[i][1]*tlist[i][1]/1.e6
	x2[i] = x[i]
	y2[i] = 2.*x2[i]

plt.plot(x,y,'rx')
plt.plot(x2,y2,'b')
plt.xlabel('Genomic distance in Mbp')
plt.ylabel('Mean square distance ($\mu m^2$)')
plt.axis([0., 200., 0., 40.])
# plt.savefig("dod.png")
plt.show()
