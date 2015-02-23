
#include <iostream>
#include <cmath>
#include <utility>
#include "replicator.hh"

using std::cerr;
using std::endl;
using std::fabs;
using std::pair;

bool firing_test()
{
	bool all_ok = true;
	
	replicator repli;
	
	parameter_set params;
	params.strand_length = 1.e7;
	
	repli.set_parameters(params);
	chromatin_type t1;
	t1.fork_speed = 28.;
	repli.chromatins.push_back(t1);
	
	boundary b1(0.,true, 0, 0);
	boundary b2(repli.parameters.strand_length, true, 0, 0);
	repli.boundaries.insert(b1);
	repli.boundaries.insert(b2);
	
	origin o1(0.2*repli.parameters.strand_length, 0);
	
	red_black_tree<origin>::iterator p1 = repli.origins.insert(o1);
	
	repli.cascade_counter = 23;
	pair<red_black_tree<origin>::iterator, unsigned int> res = repli.get_firing_origin();
	red_black_tree<origin>::iterator p2 = res.first;
	
	if (not (p1 == p2))
	{
		cerr << "Origin selector did not return the only origin available in the spontaneous firing case." << endl;
		all_ok = false;
	}
	
	repli.fire_origin(p2, res.second);
	
	if (repli.origins.get_n_elements() != 0)
	{
		cerr << "Fired origin was not removed." << endl;
		all_ok = false;
	}
	
	origin o2(0.2*repli.parameters.strand_length + 1.e4, 0);
	
	red_black_tree<origin>::iterator p3 = repli.origins.insert(o2);
	pair<red_black_tree<origin>::iterator, unsigned int> res2 = repli.get_firing_origin();
	
	red_black_tree<origin>::iterator p4 = res2.first;
	if (p3 != p4)
	{
		cerr << "Origin selector did not return the only origin available in the induced firing case." << endl;
		all_ok = false;
	}
	
	repli.fire_origin(p4, res.second);
	
	if (repli.origins.get_n_elements() != 0)
	{
		cerr << "Fired origin was not removed." << endl;
		all_ok = false;
	}
	
	if (not (res2.second == 23))
	{
		cerr << "Cascade index was not propagated correctly." << endl;
		all_ok = false;
	}
	
	return all_ok;
}

int main(int argc, char* argv[])
{
	return firing_test();
}
