#!/usr/bin/env python2
# -*- coding: utf-8 -*-

n_repititions = 10
n_entries = 20

infile = open("results/distcheck4/res_distances_9999.txt")
filelines = []
for line in infile:
	filelines.append(line)
infile.close()

valists = []
valists2 = []
for k in range(0,n_repititions):
	valists.append([])
	valists2.append([])

for j in range(0,n_entries):
	for k in range(0,n_repititions):
		valists[k].append(filelines[j*n_repititions+k])
		valists2[k].append(filelines[n_entries*n_repititions+j*n_repititions+k])
		

outfile = open("results/distcheck4/fixed_res_distances_9999.txt","w")
for k in range(0,n_repititions):
	for j in range(0,n_entries):
		outfile.write(valists[k][j])
for k in range(0,n_repititions):
	for j in range(0,n_entries):
		outfile.write(valists2[k][j])
outfile.close()