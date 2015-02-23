#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re,sys,os,subprocess
from test_helpers import resprint, RED, GREEN, CEND


fulldir = os.path.abspath(os.path.dirname(__file__))
makedir = os.path.dirname(fulldir)
tmdname = "st_test"

openstring = "rm -rf " + fulldir + "/"+ tmdname
p = os.popen(openstring)

print "Creating simtools test install... ",
os.chdir(makedir)

confstrings = ['./configure', '--disable-optimization', '--enable-debug' ,'--libdir=' + fulldir + "/" + tmdname, "--includedir=" + fulldir + "/" + tmdname]
makestrings = ['make', 'install']
pres = subprocess.Popen(confstrings, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
pres.wait()

errll = []
for line in pres.stderr:
	errll.append(line)

pres2 = subprocess.Popen(makestrings, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
pres2.wait()

for line in pres2.stderr:
	errll.append(line)

errll = []
for line in pres.stderr:
	errll.append(line)

resprint([not pres.returncode,errll])

os.chdir(fulldir)

import tree_test as ttt

print "Building tree tests... ",
if resprint(ttt.build(True)):
	failed = True

#print "Binary tree test...",
#if resprint(ttt.binary_test()):
	#failed = True

#print "Growing binary tree test...",
#if resprint(ttt.growing_test()):
	#failed = True

#print "Heap test...",
#if resprint(ttt.heap_test()):
	#failed = True

print "Red-black iterator test...",
if resprint(ttt.red_black_iterator_test()):
	failed = True