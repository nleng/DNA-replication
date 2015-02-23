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
#include <cmath>
#include <list>
#include <tr1/random>
#include <boost/generator_iterator.hpp>

#include "st_test/simtools/growing_binary_search.hh"

using namespace std;
using namespace simtools;

int main(int argc, char* argv[])
{
	long unsigned worked = 1;
	long unsigned failed = 0;
	
	bool all_ok = true;
	
	growing_binary_search<double> gbs1;
	
	long unsigned res1 = gbs1.insert(5.);
	
	if ((not check_growing_internals(gbs1, 5., 0, 1, 1)) or (res1 != 0))
	{
		cerr << "hn_growing_binary_search shows inconsistency in first check." << endl;
		all_ok = false;
	}
	long unsigned res2 = gbs1.insert(3.);
	
	if ((not check_growing_internals(gbs1, 8., 0, 2, 2)) or (res2 != 1))
	{
		cerr << "hn_growing_binary_search shows inconsistency in second check." << endl;
		all_ok = false;
	}
	
	long unsigned res3 = gbs1.insert(22.);
	
	if ((not check_growing_internals(gbs1, 30., 0, 4, 3)) or (res3 != 2))
	{
		cerr << "hn_growing_binary_search shows inconsistency in third check." << endl;
		all_ok = false;
	}
	
	gbs1.remove(1);
	
	if (not check_growing_internals(gbs1, 27., 0, 4, 3))
	{
		cerr << "hn_growing_binary_search shows inconsistency in fourth check." << endl;
		all_ok = false;
	}
	
	gbs1.remove(0);
	
	if (not check_growing_internals(gbs1, 22., 0, 4, 3))
	{
		cerr << "hn_growing_binary_search shows inconsistency in fifth check." << endl;
		all_ok = false;
	}
	
	long unsigned res4 = gbs1.insert(1.);
	
	if ((not check_growing_internals(gbs1, 23., 0, 4, 3)) or (res4 != 3))
	{
// 		cout << "r3 " << res4 << endl;
		cerr << "hn_growing_binary_search shows inconsistency in sixth check." << endl;
		all_ok = false;
	}
	
	long unsigned res5 = gbs1.insert(1.);
	
	if ((not check_growing_internals(gbs1, 24., 2, 4, 3)) or (res5 != 4))
	{
		cerr << "hn_growing_binary_search shows inconsistency in seventh check." << endl;
		all_ok = false;
	}
	
	long unsigned res6 = gbs1.insert(4.);
	
	if ((not check_growing_internals(gbs1, 28., 2, 4, 3)) or (res6 != 5))
	{
		cerr << "hn_growing_binary_search shows inconsistency in eighth check." << endl;
		all_ok = false;
	}
	
	long unsigned res7 = gbs1.insert(4.);
	
	if ((not check_growing_internals(gbs1, 32., 2, 8, 4)) or (res7 != 6))
	{
		cerr << "hn_growing_binary_search shows inconsistency in ninth check." << endl;
		all_ok = false;
	}
	
	gbs1.remove(5);
	
	if (not check_growing_internals(gbs1, 28., 2, 8, 4))
	{
		cerr << "hn_growing_binary_search shows inconsistency in tenth check." << endl;
		all_ok = false;
	}
	
	long unsigned ind1 = gbs1.get_index(24.1);
	if (ind1 != 6)
	{
		cerr << "First index retrieval test result wrong. Is " << ind1 << ", but should be " << 6 << "." << endl;
		all_ok = false;
	}
	
	long unsigned ind2 = gbs1.get_index(23.9);
	if (ind2 != 4)
	{
		cerr << "Second index retrieval test result wrong. Is " << ind2 << ", but should be " << 4 << "." << endl;
		all_ok = false;
	}
	
	long unsigned ind3 = gbs1.get_index(24.0);
	if (ind3 != 4)
	{
		cerr << "Third index retrieval test result wrong. Is " << ind3 << ", but should be " << 4 << "." << endl;
		all_ok = false;
	}
	
	// Do some fuzzing.
	std::tr1::mt19937 engine(static_cast<std::tr1::mt19937::result_type>(23));
	std::tr1::uniform_real<double> ranval(0., 1.);
	std::tr1::variate_generator<std::tr1::mt19937, std::tr1::uniform_real<double> > VG1(engine,ranval);
	boost::generator_iterator<std::tr1::variate_generator<std::tr1::mt19937, std::tr1::uniform_real<double> > > trrs1(&VG1);
	
	size_t n_test_entries = 1000;
	size_t n_remove = 500;
	size_t n_test = 1000;
	std::tr1::uniform_int<int> rancon(0, n_test_entries-1);
	
	for(size_t i = 0; i < n_test_entries; i++)
	{
		double tval = *trrs1++;
		gbs1.insert(tval);
	}
	for(size_t j = 0; j < n_remove; j++)
	{
		long unsigned tval = rancon(engine) + 2;
		gbs1.remove(tval);
	}
	activity_consistency(gbs1);
	
	double max_prop = gbs1.get_propensity();
	for(size_t k = 0; k < n_test; k++)
	{
		double tval = *trrs1++;
		tval *= max_prop;
		try
		{
			gbs1.get_index(tval);
		}
		catch(...)
		{
			cerr << "Caught something in iteration " << k << " of growing tree fuzzing test." << endl;
		}
	}
	for(size_t i = 0; i < n_test_entries; i++)
	{
		double tval = *trrs1++;
		gbs1.insert(tval);
	}
	activity_consistency(gbs1);
	
	if (not all_ok)
		return failed;
	return worked;
}