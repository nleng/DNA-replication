#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#infile = open("results/res_distances_5000.txt")
#infile = open("results/distcheck4/res_distances_9999.txt")
#infile = open("results/distcheck5/res_distances_3150.txt")
#infile = open("results/distcheck6/res_distances_2960.txt")
#infile = open("results/distcheck7/res_distances_2950.txt")
infile = open("/imports/bott.work/loeb/dna_metropolis/distckeck_103/res_distances_99.txt")
#infile = open("/imports/bott.work/loeb/dna_metropolis/distckeck_104/res_distances_99.txt")
#infile = open("/imports/bott.work/loeb/dna_metropolis/distckeck_105/res_distances_99.txt")
#infile = open("/imports/bott.work/loeb/dna_metropolis/distckeck_106/res_distances_99.txt")
infile = open("/imports/bott.work/loeb/dna_metropolis/distckeck_108/res_distances_99.txt")

n_entries = 60
n_repititions = 2*1200/n_entries

distlist = []
xlist = []
distlist2 = []
xlist2 = []
cntr = 0
for line in infile:
	zz = line.replace('\n','').split('\t')
	
	if cntr < n_entries*n_repititions: 
		xlist.append(float(zz[1])-float(zz[0]))
		distlist.append(float(zz[2]))
	else:
		xlist2.append(float(zz[1])-float(zz[0]))
		distlist2.append(float(zz[2]))
	cntr += 1

x = np.zeros(n_entries)
y = np.zeros(n_entries)
x2 = np.zeros(n_entries)
y2 = np.zeros(n_entries)

maxv = 0.
xmaxv = 0.
for i in range(0,n_entries):
	for k in range(0,n_repititions):
	#for k in range(0,1):
		if k == 0:
			x[i] = xlist[i]/1.e6
		y[i] += (distlist[i+k*n_entries]/1.e3)*(distlist[i+k*n_entries]/1.e3)/n_repititions
		#y[i] += (distlist[i+k*n_entries]/1.e3)*(distlist[i+k*n_entries]/1.e3)
		if y[i] > maxv:
			maxv = y[i]
		if x[i] > xmaxv:
			xmaxv = x[i]
		
		
for i in range(0,n_entries):
	for k in range(0,n_repititions):
	#for k in range(0,1):
		if k == 0:
			x2[i] = xlist2[i]/1.e6
		y2[i] += (distlist2[i+k*n_entries]/1.e3)*(distlist2[i+k*n_entries]/1.e3)/n_repititions
		#y2[i] += (distlist2[i+k*n_entries]/1.e3)*(distlist2[i+k*n_entries]/1.e3)
		if i == 30:
                        print(distlist2[i+k*n_entries]/1.e3)

		if y2[i] > maxv:
			maxv = y2[i]
		if x2[i] > xmaxv:
			xmaxv = x2[i]

fg1 = plt.figure()

plt.plot(x,y,'bx')
plt.plot(x,y,'g',label="euchromatin")
plt.plot(x2,y2,'bx')
plt.plot(x2,y2,'r',label="heterochromatin")
plt.xlabel('Genomic distance (Mbp)')
plt.ylabel('Mean sq. displacement ($\mu m^2$)')
plt.legend()
plt.axis([0, xmaxv, 0, maxv])
plt.show()
