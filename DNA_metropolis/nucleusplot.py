#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys

from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

import FTGL

from copy import copy
from math import sin,cos, pi
from libreplihelpers import get_3d_points

#if os.getenv('USER') == "spock":
	#vertfile = "/home/spock/programming/replication/src/results/res_fork_0.0333333_400000_67000_5e+06_500000_0.txt"
#elif os.getenv('USER') == "garak":
	##vertfile = "results/paper_res_fork_0.txt"
	##beadfilename = "results/fullcalc14.txt"
	##beadfilename = "results/fullcalc15.txt"
	#vertfile = "/home/garak/phd/programs/replication/src/results/paper_hela_realistic_forks/res_fork_6000_280000_67000_28_500000_0_run.txt"
	#beadfilename = "results/paper_random/res_res_999.txt"
	##beadfilename = "results/loops1.txt"
#elif os.getenv('USER') == "loeb":
##	vertfile = "/home/loeb/phd/programs/replication/src/results/res_fork_0.0333333_400000_67000_5e+06_500000_0.txt"
	#vertfile = "res_fork_6000_280000_67000_28_500000_0_run.txt"

	##beadfilename = "paper_res.txt"
	#beadfilename = "/imports/bott.work/loeb/dna_metropolis/paper_alterna78/res_res_390.txt"
	##beadfilename = "/imports/hus.work/loeb/dna_metropolis/loops3/res_res_9999.txt"

vertfile = '/home/corni/aaReplication/replicationSim/results/paper_hela_realistic_forks/res_fork_6000_280000_67000_28_500000_0_run.txt'
beadfilename = '/home/corni/aaReplication/dnaSim/outFile/res_res_12.txt'
#with_color = True
#onlydots = True
radius = 7500.
#framenumber=350

paramslist = [	["D", 90, False, False],
		["E", 180, False, False],
		["F", 351, False, False],
		["G", 90, True, False],
		["H", 180, True, False],
		["I", 351, True, True]
	]


number = 5


letter = paramslist[number][0]
onlydots = True
#radius = 7500.
#framenumber=351
framenumber= 100	# paramslist[number][1]
#framenumber=180
usecolor = paramslist[number][2]
plotlegend = paramslist[number][3]


def initFun():
	global font
	global font2
	glClearColor(0.0, 0.0, 0.0,0.0)
	glColor3f(1.0, 1.0, 1.0)
	glPointSize(4.0)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	gluOrtho2D(-1.1*radius,1.1*radius,-1.1*radius,1.1*radius)
	
	fpath = os.path.abspath("../../fonts/Arial.ttf")
	#print(fpath)
	font = FTGL.PolygonFont(fpath)
	font.FaceSize(24)
	font.UseDisplayList(True)
	
	fpath2 = os.path.abspath("../../fonts/Arial_Bold.ttf")
	font2 = FTGL.PolygonFont(fpath2)
	font2.FaceSize(24)
	font2.UseDisplayList(True)
	

def printtext(x, y, String):
	glMatrixMode(GL_PROJECTION)
	glPushMatrix()
	glLoadIdentity()
	glOrtho(0, 640, 0, 480, -1.0, 1.0)
	glMatrixMode(GL_MODELVIEW)
	glPushMatrix()
	glLoadIdentity()
	glPushAttrib(GL_DEPTH_TEST)
	glDisable(GL_DEPTH_TEST)
	glRasterPos2i(x,y)
	for i in range(0, len(String)):
		print(String[i])
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ord(String[i]))
	
	glPopAttrib()
	glMatrixMode(GL_PROJECTION)
	glPopMatrix()
	glMatrixMode(GL_MODELVIEW)
	glPopMatrix()

def printtext2(x, y, String, fnt):
	glMatrixMode(GL_PROJECTION)
	glPushMatrix()
	glLoadIdentity()
	glOrtho(0, 640, 0, 480, -1.0, 1.0)
	glMatrixMode(GL_MODELVIEW)
	glPushMatrix()
	glLoadIdentity()
	glPushAttrib(GL_DEPTH_TEST)
	glDisable(GL_DEPTH_TEST)
	glRasterPos2i(x,y)
	glColor3f(1.0, 1.0, 1.0)
	
	glTranslatef(x, y, 0.0)
	##glTranslatef(x, y, 0.0)
	fnt.Render(String)
	
	#for i in range(0, len(String)):
		#print(String[i])
		#glutBitmapCharacter(GLUT_BITMAP_9_BY_15, ord(String[i]))
	
	glPopAttrib()
	glMatrixMode(GL_PROJECTION)
	glPopMatrix()
	glMatrixMode(GL_MODELVIEW)
	glPopMatrix()

	
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

def displayFun():
	global rvertlist
	global gvertlist
	global font
	global font2
	
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
	
	hours = int(float(framenumber)/36.)
	minutes = int((float(framenumber)-float(hours)*36.)/0.6)
	
	# ATTENTION! Although you can not see it, these spaces are unicode thin spaces.
	
	toshowstr = "t="+str(hours) +"h "
	if minutes != 0:
		toshowstr += str(minutes) + "min"
	
	printtext2(10,450, letter, font2)
	printtext2(40,440, toshowstr, font)
	
	if plotlegend:
		glColor(255,255,255)
		glLineWidth(2.)
		glBegin(GL_LINES)
		glVertex2f(5300,-6800)
		glVertex2f(7300, -6800)
		
		glEnd()
		glFlush()
		
		printtext2(540, 20, "2μm", font)
		
	
	glutSwapBuffers()
	#font.Render("foo")
	#glColor3f(.0, 1.0, 1.0)
	#glColor3f(1.0, 1.0, 1.0)
	#x = 0
	#y = 0
        
	##glPushMatrix()
	##glTranslatef(x, y, 0.0)
	#font.Render("foo")
	
	#glFlush()
	#glPopMatrix()
	#glutSwapBuffers()
	
        #for j in range(0, 4):
            #y = 275.0 - i * 120.0 - j * yild
            #if i >= 3:
                #glRasterPos(x, y)
                #font.Render(string[j])
            #elif i == 2:
                #glEnable(GL_TEXTURE_2D)
                #glEnable(GL_BLEND)
                #glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
                #glPushMatrix()
                #glTranslatef(x, y, 0.0)
                #font.Render(string[j])
                #glPopMatrix()
                #glDisable(GL_TEXTURE_2D)
                #glDisable(GL_BLEND)
            #else:
	

if __name__ == '__main__':
	global rvertlist
	global gvertlist
	global font
	
	if not (os.path.exists(vertfile)):
		print("Fork file not found.")
		sys.exit()
	
	[gvertlist, rvertlist] = get_vertices_from_file(vertfile, framenumber)
      	
	#if not with_color:
		#gvertlist += rvertlist
		#rvertlist = []
	print("N_eu:", len(gvertlist), "N_het:", len(rvertlist))
	
	if not usecolor:
		gvertlist += rvertlist
		rvertlist = []
	glutInit()
	glutInitWindowSize(640,480)
	glutCreateWindow("Drawdots")
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB)
	glutDisplayFunc(displayFun)
	initFun()
	glutMainLoop()
