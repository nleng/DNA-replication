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

template<typename data_t> bool checkrbt(const red_black_tree<data_t> &rbt)
{
	bool all_ok = true;
	
	pair<bool, string>  res1 = rbt.check_sortedness();
	if (not res1.first)
	{
		cerr << res1.second << endl;
		all_ok = false;
	}
	
	pair<bool, string> res2 = rbt.check_family_sizes();
	if (not res2.first)
	{
		cerr << res2.second << endl;
		all_ok = false;
	}
	
	pair<bool, string> res3 = rbt.check_root_blackness();
	if (not res3.first)
	{
		cerr << res3.second << endl;
		all_ok = false;
	}
	
	pair<bool, string> res4 = rbt.check_color_correctness();
	if (not res4.first)
	{
		cerr << res4.second << endl;
		all_ok = false;
	}
	
	pair<bool, string> res5 = rbt.check_black_count();
	if (not res5.first)
	{
		cerr << res5.second << endl;
		all_ok = false;
	}
	
	pair<bool, string> res6 = rbt.check_backlink_consistency();
	if (not res6.first)
	{
		cerr << res6.second << endl;
		all_ok = false;
	}
	return all_ok;
}

int main(int argc, char* argv[])
{
	bool all_ok = true;
	
	std::tr1::mt19937 engine(static_cast<std::tr1::mt19937::result_type>(23));
	
	std::tr1::uniform_real<double> valdist(0., 100.);
	std::tr1::variate_generator<std::tr1::mt19937, std::tr1::uniform_real<double> > VG1(engine,valdist);
	boost::generator_iterator<std::tr1::variate_generator<std::tr1::mt19937, std::tr1::uniform_real<double> > > trrs1(&VG1);
	
	red_black_tree<double> trbt;
	trbt.insert(3.);
	trbt.insert(5.);
	trbt.insert(1.);
	trbt.insert(12.);
	trbt.insert(77.);
	
	if (not checkrbt<double>(trbt))
		all_ok = false;
	
	for (size_t i = 0; i < 1000; i++)
	{
		trbt.insert(*trrs1++);
	}
	if (not checkrbt<double>(trbt))
		all_ok = false;
	
	size_t n_entries = trbt.get_n_elements();
	
	for (size_t i = 0; i < n_entries; i++)
	{
// 		cout << i << endl;
		red_black_tree<double>::iterator p = trbt.find_smallest();
		trbt.remove(p);
		if (not checkrbt<double>(trbt))
		{
			all_ok = false;
			break;
		}
	}	
	
	if (not checkrbt<double>(trbt))
		all_ok = false;
	
	for (size_t i = 0; i < 1000; i++)
	{
		trbt.insert(*trrs1++);
	}
	
	for (size_t i = 0; i < 1000; i++)
	{
		uniform_int_t my_uni(0, trbt.get_n_elements()-1);
		uni_int_gen ug(engine, my_uni);
		boost::generator_iterator<uni_int_gen> ugi(&ug);
		unsigned long tmp = *ugi++;
		red_black_tree<double>::iterator p = trbt.find_smallest(tmp);
		trbt.remove(p);
		
		if (not checkrbt<double>(trbt))
		{
			all_ok = false;
			break;
		}
	}
	
	if (not (trbt.get_n_elements() == 0))
	{
		cerr << "Red black tree content assert failed." << endl;
		all_ok = false;
	}
	
	for (size_t i = 0; i < 1000; i++)
	{
		trbt.insert(*trrs1++);
	}
	
	double sval = *(trbt.find_smallest());
	double bval = *(trbt.find_biggest());
	
	red_black_tree<double>::iterator low = trbt.find_smaller_than(sval/2.);
	red_black_tree<double>::iterator high = trbt.find_bigger_than(bval*2.);
	
	if (not((low == NULL) or (high == NULL)))
	{
		cerr << "Red black tree returned iterator for out of bound request." << endl;
		all_ok = false;
	}
	
	red_black_tree<double>::iterator low2 = trbt.find_smaller_than((sval+bval)/2.);
	red_black_tree<double>::iterator high2 = trbt.find_bigger_than((sval+bval)/2.);
	
	if ((low2 == NULL) or (high2 == NULL))
	{
		cerr << "Red black tree failed to return iterator for in bound request." << endl;
		all_ok = false;
	}
	
	if ((not (*low2 < (sval+bval)/2.)) or (not (*high2 > (sval+bval)/2.)))
	{
		cerr << "Red black tree failed to value lower or higher than requested value." << endl;
		all_ok = false;
	}
	
	return all_ok;
}