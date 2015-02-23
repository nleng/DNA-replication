// simtools - A C++ toolset for stochastic network dynamics models.
// 
// Copyright © 2010-2012 Daniel Löb <daniel@zombiepiratesfromspace.eu>
// 
// simtools is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// simtools is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with simtools.  If not, see <http://www.gnu.org/licenses/>.


#include <iostream>
#include <tr1/random>
#include <boost/generator_iterator.hpp>
#include "st_test/simtools/red_black_tree.hh"

using std::cout;
using std::endl;
using std::cerr;
using std::endl;
using namespace simtools;

typedef std::tr1::uniform_int<unsigned long> uniform_int_t;
typedef std::tr1::variate_generator<std::tr1::mt19937, uniform_int_t> uni_int_gen;

int main(int argc, char* argv[])
{
	bool all_ok = true;
	
	red_black_tree<double> trbt;
	trbt.insert(40.5);
	trbt.insert(6.);
	trbt.insert(0.5);
	trbt.insert(22.);
	trbt.insert(73.);
	
	red_black_tree<double>::iterator p = trbt.find_smallest();
	if(not (*p == 0.5))
	{
		cerr << "Smallest position value is wrong." << endl;
		all_ok = false;
	}
	
	*p = 0.7;
	
	red_black_tree<double>::iterator p2 = trbt.find_smallest();
	if(not (*p2 == 0.7))
	{
		cerr << "Changing of smallest value failed." << endl;
		all_ok = false;
	}
	
	red_black_tree<pair<double,size_t> > xxx;
	xxx.insert(pair<double,size_t>(0.3,12));
	xxx.insert(pair<double,size_t>(23.,55));
	xxx.insert(pair<double,size_t>(12.,3));
	
	red_black_tree<pair<double,size_t> >::iterator pit = xxx.find_smallest();
	
	if(not (pit->first == 0.3))
	{
		cerr << "Smallest position value is wrong in pair test." << endl;
		all_ok = false;
	}
	
	pit->first = 0.8;
	
	if(not (pit->first == 0.8))
	{
		cerr << "Changing of smalles pair value failed." << endl;
		all_ok = false;
	}
	
	return all_ok;
}