#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os, sys

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

from copy import copy
from math import sin,cos, pi
from libreplihelpers import get_3d_points

if os.getenv('USER') == "spock":
	vertfile = "/home/spock/programming/replication/src/results/res_fork_0.0333333_400000_67000_5e+06_500000_0.txt"
elif os.getenv('USER') == "garak":
	vertfile = "res_fork_0.0333333_400000_67000_5e+06_500000_0.txt"
	#beadfilename = "dna_metropolis/fixed_twozone/res_res_611.txt"
	#beadfilename = "dna_metropolis/fixed_benchmark/res_res_5850.txt"
	beadfilename = "results/res_res_5670.txt"
	#beadfilename = "dna_metropolis/fixed_benchmark2/res_res_5860.txt"
elif os.getenv('USER') == "loeb":
	vertfile = "/home/loeb/phd/programs/replication/src/results/res_fork_0.0333333_400000_67000_5e+06_500000_0.txt"
#	beadfilename = "/imports/lydia.work/loeb/dna_metropolis/ranlo8_8/res_res_9999.txt"
#	beadfilename = "/imports/lydia.work/loeb/dna_metropolis/benchmark4/res_res_1962.txt"
#	beadfilename = "/imports/lydia.work/loeb/dna_metropolis/benchmark6_c4/res_res_2520.txt"
	#beadfilename = "/imports/bott.work/loeb/dna_metropolis/twozone3/res_res_8800.txt"
	beadfilename = "/imports/debora.work/loeb/dna_metropolis/twozone_7/res_res_2470.txt"

onlydots = True
radius = 5000.
framenumber=100

def initFun():
	glClearColor(0.0, 0.0, 0.0,0.0)
	glColor3f(1.0, 1.0, 1.0)
	glPointSize(4.0)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	gluOrtho2D(-1.1*radius,1.1*radius,-1.1*radius,1.1*radius)


def get_vertices_from_file(fname, toshow):
	infile = open(fname)
	
	counter = 0
	greendatalist = []
	reddatalist = []
	for line in infile:
		if counter == toshow:
			zz = line.replace('\n','').split('\t')
			for i in range(1,len(zz)):
				#if float(zz[i]) < 225949719.0:
				if float(zz[i]) < 60.e6:
					greendatalist.append(float(zz[i]))
				elif float(zz[i]) < 120.e6:
					reddatalist.append(float(zz[i]))
			greendatalist.sort()
			reddatalist.sort()
		counter += 1
		
	gtmplist = get_3d_points(greendatalist, beadfilename)
	rtmplist = get_3d_points(reddatalist, beadfilename)
	
	print("done")
	return [gtmplist, rtmplist]

def displayFun():
	global rvertlist
	global gvertlist
	
	glClear(GL_COLOR_BUFFER_BIT)
	# Draw the Cell itself
	tvlist = []
	#pi = 3.14159
	maxp = 200
	for i in range(0, maxp):
		tvlist.append([radius*sin(pi*i*2./maxp),radius*cos(pi*i*2./maxp)])
	
	glColor(255,255,255)
	glLineWidth(2.)
	glBegin(GL_LINES)
	for i in range(0,len(tvlist)-1):
		glVertex2f(tvlist[i][0], tvlist[i][1])
		glVertex2f(tvlist[i+1][0], tvlist[i+1][1])
	glVertex2f(tvlist[len(tvlist)-1][0], tvlist[len(tvlist)-1][1])
	glVertex2f(tvlist[0][0], tvlist[0][1])
	glEnd()
	glFlush()
	# Draw the replication positions
	if onlydots:
		glColor(0,255,0)
		glPointSize(2.)
		glBegin(GL_POINTS)
		for i in range(0, len(gvertlist)-1):
			#print vertlist[i]
			#glVertex3f(vertlist[i][0], vertlist[i][1], vertlist[i][2])
			#glVertex3f(vertlist[i+1][0], vertlist[i+1][1], vertlist[i+1][2])
			
			glVertex2f(gvertlist[i][0], gvertlist[i][1])
			#glVertex2f(vertlist[i+1][0], vertlist[i+1][1])
		glEnd()
		glFlush()
		
		
		glColor(255,0,0)
		glPointSize(2.)
		glBegin(GL_POINTS)
		for i in range(0, len(rvertlist)-1):
			#print vertlist[i]
			#glVertex3f(vertlist[i][0], vertlist[i][1], vertlist[i][2])
			#glVertex3f(vertlist[i+1][0], vertlist[i+1][1], vertlist[i+1][2])
			
			glVertex2f(rvertlist[i][0], rvertlist[i][1])
			#glVertex2f(vertlist[i+1][0], vertlist[i+1][1])
		glEnd()
		glFlush()
	else:
		#vertlist = get_vertex_list(1000000,100.)
		glClear(GL_COLOR_BUFFER_BIT)
		glBegin(GL_LINES)
		for i in range(0, len(vertlist)-1):
			#print vertlist[i]
			#glVertex3f(vertlist[i][0], vertlist[i][1], vertlist[i][2])
			#glVertex3f(vertlist[i+1][0], vertlist[i+1][1], vertlist[i+1][2])
			
			glVertex2f(vertlist[i][0], vertlist[i][1])
			glVertex2f(vertlist[i+1][0], vertlist[i+1][1])
		glEnd()
		glFlush()

if __name__ == '__main__':
	global rvertlist
	global gvertlist 
	[gvertlist,rvertlist] = get_vertices_from_file(vertfile, framenumber)
	glutInit()
	glutInitWindowSize(640,480)
	glutCreateWindow("Drawdots")
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB)
	glutDisplayFunc(displayFun)
	initFun()
	glutMainLoop()
