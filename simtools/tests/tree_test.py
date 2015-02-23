#!/usr/bin/env python
# -*- coding: utf-8 -*-
# simtools - A C++ toolset for stochastic network dynamics models.
# 
# Copyright © 2010-2012 Daniel Löb <daniel@zombiepiratesfromspace.eu>
# 
# simtools is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# simtools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with simtools.  If not, see <http://www.gnu.org/licenses/>.

import os
import subprocess
from test_helpers import runtest

def build(debug=False):
	errll = []
	apath = os.path.abspath(os.path.dirname(__file__))
	curdir = apath+"/st_test"
	
	failed = False
	
	filename = 'binary_test'
	
	toremove = apath + filename
	if(os.path.exists(toremove)):
		os.remove(toremove)
	
	if debug:
		compilecommands = ['g++', '-Wall', '-ggdb3', '-std=c++11', '-L' +curdir, '-lsimtools', 'bs_test.cpp', '-o', filename]
	else:
		compilecommands = ['g++', '-Wall', '-O3', '-std=c++11', '-L' +curdir, '-lsimtools', 'bs_test.cpp', '-o', filename]
	
	pres = subprocess.Popen(compilecommands, stderr=subprocess.PIPE)
	pres.wait()
	
	for line in pres.stderr:
		errll.append(line)
	
	filename = 'growing_test'
	
	toremove = apath + filename
	if(os.path.exists(toremove)):
		os.remove(toremove)
		
	if debug:
		compilecommands = ['g++', '-Wall', '-ggdb3', '-std=c++11', '-L' +curdir, '-lsimtools', 'growing_test.cpp', '-o', filename]
	else:
		compilecommands = ['g++', '-Wall', '-O3', '-std=c++11', '-L' +curdir, '-lsimtools', 'growing_test.cpp', '-o', filename]
	
	pres2 = subprocess.Popen(compilecommands, stderr=subprocess.PIPE)
	pres2.wait()
	
	for line in pres2.stderr:
		errll.append(line)
	
	filename = 'heap_test'
	
	toremove = apath + filename
	if(os.path.exists(toremove)):
		os.remove(toremove)
		
	if debug:
		compilecommands = ['g++', '-Wall', '-ggdb3', '-std=c++11', '-L' +curdir, '-lsimtools', 'heap_test.cpp', '-o', filename]
	else:
		compilecommands = ['g++', '-Wall', '-O3', '-std=c++11', '-L' +curdir, '-lsimtools', 'heap_test.cpp', '-o', filename]
	
	pres3 = subprocess.Popen(compilecommands, stderr=subprocess.PIPE)
	pres3.wait()
	
	for line in pres3.stderr:
		errll.append(line)
	
	filename = 'rb_test'
	
	toremove = apath + filename
	if(os.path.exists(toremove)):
		os.remove(toremove)
		
	if debug:
		compilecommands = ['g++', '-std=c++11', '-Wall', '-ggdb3', '-L' +curdir, '-lsimtools', 'rb_test.cpp', '-o', filename]
	else:
		compilecommands = ['g++', '-std=c++11', '-Wall', '-O3', '-L' +curdir, '-lsimtools', 'rb_test.cpp', '-o', filename]
	
	pres4 = subprocess.Popen(compilecommands, stderr=subprocess.PIPE)
	pres4.wait()
	
	for line in pres4.stderr:
		errll.append(line)
	
	filename = 'rb_iterator_test'
	
	toremove = apath + filename
	if(os.path.exists(toremove)):
		os.remove(toremove)
		
	if debug:
		compilecommands = ['g++', '-std=c++11', '-Wall', '-ggdb3', '-L' +curdir, '-lsimtools', 'rb_iterator_test.cpp', '-o', filename]
	else:
		compilecommands = ['g++', '-std=c++11', '-Wall', '-O3', '-L' +curdir, '-lsimtools', 'rb_iterator_test.cpp', '-o', filename]
	
	pres5 = subprocess.Popen(compilecommands, stderr=subprocess.PIPE)
	pres5.wait()
	
	for line in pres5.stderr:
		errll.append(line)
	
	if len(errll) > 0:
		failed = True
	
	return [not failed, errll]

def growing_test():
	return runtest('./growing_test')

def binary_test():
	return runtest('./binary_test')

def heap_test():
	return runtest('./heap_test')

def red_black_test():
	return runtest('./rb_test')

def red_black_iterator_test():
	return runtest('./rb_iterator_test')
