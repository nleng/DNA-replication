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

#include "st_test/simtools/binary_heap.hh"
#include <iostream>
#include <cmath>
#include <tr1/random>
#include <boost/generator_iterator.hpp>

using std::cout;
using std::endl;
using std::cerr;
using std::abs;
using namespace simtools;

int main(int argc, char* argv[])
{
	long unsigned worked = 1;
	long unsigned failed = 0;
	
	bool all_ok = true;
	
	binary_heap<double> mybin;
	
	mybin.insert_value(4.,0);
	mybin.insert_value(3.,1);
	mybin.insert_value(3.5,2);
	
	if ((mybin.get_smallest()).second != 1)
	{
		cerr << "Wrong smallest value in heap returned. Has index " << (mybin.get_smallest()).second << " but should be " << 1 << "."<< endl;
		all_ok = false;
	}
	
	mybin.insert_value(5.,3);
	mybin.insert_value(1.,4);
	mybin.insert_value(3.,5);
	
	if ((mybin.get_smallest()).second != 4)
	{
		cerr << "Wrong smallest value in heap returned. Has index " << (mybin.get_smallest()).second << " but should be " << 4 << "."<< endl;
		all_ok = false;
	}
	
	if (abs(mybin.get_value(5) - 3.) > 1e-10)
	{
		cerr << "Wrong value returned in first get value. Is " << mybin.get_value(5) << " but should be " << 3. << "." << endl;
		all_ok = false;
	}
	
	long unsigned sval = (mybin.remove_smallest()).second;
	
	if (sval != 4)
	{
		cerr << "Wrong smallest value removed from heap. Has index " << sval << " but should be " << 4 << "." << endl;
		all_ok = false;
	}
	
	if ((mybin.get_smallest()).second != 1)
	{
		cerr << "Wrong smallest value in heap returned. Has index " << (mybin.get_smallest()).second << " but should be " << 1 << "."<< endl;
		all_ok = false;
	}
	
	if (abs(mybin.get_value(4) + 1.) > 1e-10)
	{
		cerr << "Wrong value returned in first get value. Is " << mybin.get_value(5) << " but should be " << -1. << "." << endl;
		all_ok = false;
	}
	
	mybin.remove_smallest();
	long unsigned sval3 = (mybin.remove_smallest()).second;
	
	if (sval3 != 5)
	{
		cerr << "Wrong smallest value removed from heap. Has index " << sval3 << " but should be " << 5 << "." << endl;
		all_ok = false;
	}
	
	if (mybin.get_n_elements() != 3)
	{
		cerr << "Wrong number of elements in heap, is  " << mybin.get_n_elements() << " but should be " << 3 << "." << endl;
		all_ok = false;
	}
	
	// Now we have a situation in which there is one completely empty depth layer. After adding another element the depth should still be the same.
	mybin.insert_value(3.,6);
	
	if (mybin.get_depth() != 3)
	{
		cerr << "Wrong heap depth, is  " << mybin.get_depth() << " but should be " << 3 << "." << endl;
		all_ok = false;
	}
	
	mybin.insert_value(3.2,7);
	mybin.remove(7);
	
	// Now do some fuzzing.
	std::tr1::mt19937 engine(static_cast<std::tr1::mt19937::result_type>(23));
	std::tr1::uniform_real<double> ranval(0.1, 10.);
	std::tr1::variate_generator<std::tr1::mt19937, std::tr1::uniform_real<double> > VG1(engine,ranval);
	boost::generator_iterator<std::tr1::variate_generator<std::tr1::mt19937, std::tr1::uniform_real<double> > > trrs1(&VG1);
	
	long unsigned counter = 1000;
	for (long unsigned i = 0; i < 1; i++)
	{
		for (long unsigned j = 0; j < 1000; j++)
		{
			mybin.insert_value(*trrs1++,counter);
			counter++;
		}
		
		if (not mybin.test_consistency())
			all_ok = false;
		
		for (long unsigned j = 0; j < 1000; j++)
		{
			mybin.remove_smallest();
		}
		
		if (not mybin.test_consistency())
			all_ok = false;
	}
	
	// Test random removal.
	
	long unsigned counter2 = 1e7;
	long unsigned counter3 = 1e7;
	
	std::tr1::uniform_int<int> rancon(0, 999);
	
	for (long unsigned i = 0; i < 1; i++)
	{
		for (long unsigned j = 0; j < 1000; j++)
		{
			mybin.insert_value(*trrs1++,counter2);
			counter2++;
		}
		mybin.test_consistency();
		
		// Create a list of shuffled indentifiers.
		vector<long unsigned> rem_order(1000);
		for (long unsigned j = 0; j < 1000; j++)
		{
			rem_order[j] = counter3;
			counter3++;
		}
		
		for (long unsigned j = 0; j < 3000; j++)
		{
			size_t i1 = rancon(engine);
			size_t i2 = rancon(engine);
			long unsigned tmp = rem_order[i1];
			rem_order[i1] = rem_order[i2];
			rem_order[i2] = tmp;
		}
		
		// Remove the values again, according to this list.
		for (long unsigned j = 0; j < 1000; j++)
		{
			mybin.remove(rem_order[j]);
		}
		
		if (not mybin.test_consistency())
		{
			cerr << "Binary tree is incosistent after removal test." << endl;
			all_ok = false;
		}
	}
	
	if (not all_ok)
		return failed;
	return worked;
}