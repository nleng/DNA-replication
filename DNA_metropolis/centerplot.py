#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os, sys

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

from copy import copy
from math import sin,cos, pi

onlydots = True
radius = 5000.

def initFun():
	glClearColor(0.0, 0.0, 0.0,0.0)
	glColor3f(1.0, 1.0, 1.0)
	glPointSize(4.0)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	gluOrtho2D(-1.1*radius,1.1*radius,-1.1*radius,1.1*radius)


def get_vertices_from_file():
	#chromosomesizes = {1}
	#chromosome_frequency = {76}
	
	fname = "dna_metropolis/separate_0/res_res_999.txt"
	infile = open(fname)
	
	datalist = []
	for line in infile:
		zz = line.replace('\n','').split('\t')
		tmplist = [0., 0., 0.]
		ncnt = 0
		for i in range(0,len(zz)):
			if (i%3) == 0:
				ncnt += 1
				tmplist[0] += float(zz[i])
			
			elif (i%3) == 1:
				tmplist[1] += float(zz[i])
			
			elif (i%3) == 2:
				tmplist[2] += float(zz[i])
		
		tmplist[0] = tmplist[0]/ncnt
		tmplist[1] = tmplist[1]/ncnt
		tmplist[2] = tmplist[2]/ncnt
		
		datalist.append(tmplist)
	
	print("done")
	return datalist

def displayFun():
	global vertlist
	
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
		for i in range(0, len(vertlist)-1):
			#print vertlist[i]
			#glVertex3f(vertlist[i][0], vertlist[i][1], vertlist[i][2])
			#glVertex3f(vertlist[i+1][0], vertlist[i+1][1], vertlist[i+1][2])
			
			glVertex2f(vertlist[i][0], vertlist[i][1])
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
	global vertlist 
	vertlist = get_vertices_from_file()
	glutInit()
	glutInitWindowSize(640,480)
	glutCreateWindow("Drawdots")
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB)
	glutDisplayFunc(displayFun)
	initFun()
	glutMainLoop()
