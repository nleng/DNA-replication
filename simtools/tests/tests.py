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


import re,sys,os,subprocess
from test_helpers import resprint, RED, GREEN, CEND

fulldir = os.path.abspath(os.path.dirname(__file__))
makedir = os.path.dirname(fulldir)
tmdname = "st_test"

openstring = "rm -rf " + fulldir + "/"+ tmdname
p = os.popen(openstring)

print("Creating simtools test install...",)
os.chdir(makedir)
confstrings = ['./configure', '--libdir=' + fulldir + "/" + tmdname, "--includedir=" + fulldir + "/" + tmdname]
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

resprint([not pres.returncode,errll])

os.chdir(fulldir)

import tree_test as ttt

failed = False

print("Building tree tests... ",)
if resprint(ttt.build()):
	failed = True

print("Binary tree test...",)
if resprint(ttt.binary_test()):
	failed = True

print("Growing binary tree test...",)
if resprint(ttt.growing_test()):
	failed = True

print("Heap test...",)
if resprint(ttt.heap_test()):
	failed = True

print("Red black tree test...",)
if resprint(ttt.red_black_test()):
	failed = True

print("Red-black iterator test...",)
if resprint(ttt.red_black_iterator_test()):
	failed = True

print
if failed:
	print(RED + "There were failed Tests." + CEND)
	failed = True
else:
	print(GREEN + "All tests were successful." + CEND)
