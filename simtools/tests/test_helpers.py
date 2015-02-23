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

GREEN = '\033[92m'
RED = '\033[91m'
CEND = '\033[0m'


def resprint(sim_res):
	failed = False
	
	if not sim_res[0]:
		failed = True
		print(RED + "failed." + CEND)
		print()
		for entry in sim_res[1]:
			print('\t',entry,)
		print()
	else:
		print(GREEN + "success!" + CEND)
	
	return failed

def runtest(exec_name):
	curdir = os.path.abspath(os.path.dirname(__file__))+"/st_test"
	failed = False
	env_mapping = os.environ
	if 'LD_LIBRARY_PATH' in env_mapping:
		env_mapping['LD_LIBRARY_PATH'] = curdir + ":" + env_mapping['LD_LIBRARY_PATH']
	else:
		env_mapping['LD_LIBRARY_PATH'] = curdir
	
	execcommand = exec_name.split('\t')
	pres = subprocess.Popen(execcommand, stderr=subprocess.PIPE, stdout=subprocess.PIPE, env=env_mapping)
	pres.wait()
	
	errll = []
	for line in pres.stderr:
		failed = True
		errll.append(line)
	
	for line in pres.stdout:
		failed = True
		errll.append(line)
	
	if pres.returncode != True:
		failed = True
	
	return [not failed, errll]