
#include <iostream>
#include <cmath>
#include "replicator.hh"

using std::cerr;
using std::endl;
using std::fabs;

bool step_test()
{
	bool all_ok = true;
	
	replicator repli;
	
	parameter_set params;
	params.strand_length = 1.e7;
	
	repli.set_parameters(params);
	
	chromatin_type t1;
	t1.fork_speed = 28.;
	chromatin_type t2;
	t2.fork_speed = 50.;
	repli.chromatins.push_back(t1);
	repli.chromatins.push_back(t2);
	
	boundary b1(0.,true, 0, 0);
	boundary b2(repli.parameters.strand_length/2., false, 0, 1);
	boundary b3(repli.parameters.strand_length, true, 1, 0);
	repli.boundaries.insert(b1);
	repli.boundaries.insert(b2);
	repli.boundaries.insert(b3);
	
	origin o1(0.45*repli.parameters.strand_length, 0);
	
	red_black_tree<origin>::iterator p1 = repli.origins.insert(o1);
	
	repli.fire_origin(p1, 0);
	
	repli.step();
	
	if (repli.annihilations.get_n_elements() != 2)
	{
		cerr << "Wrong number of forks after first test step." << endl;
		all_ok = false;
	}
	
	annihilation a1 = repli.annihilations.get_value(1);
	annihilation a2 = repli.annihilations.get_value(2);
	
	if (not ((a1.at == SOLO_REMOVAL) and (not a1.left)))
	{
		cerr << "Wrong annihilation type for leftgoing fork." << endl;
		all_ok = false;
	}
	
	if (not((a2.at == SOLO_REMOVAL) and (a2.left)))
	{
		cerr << "Wrong annihilation type for boundary transition replaced fork." << endl;
		all_ok = false;
	}
	
	if (not (fabs(a2.endtime - (0.5*repli.parameters.strand_length/t2.fork_speed + repli.current_time)) < 1e-10))
	{
		cerr << "Wrong time for annihilation of replaced fork." << endl;
		all_ok = false;
	}
	
	origin o2(0.2*repli.parameters.strand_length, 0);
	repli.origins.insert(o2);
	repli.new_pairs.insert_value(repli.current_time+0.5, 0);
	
	repli.step();
	
	if (repli.annihilations.get_n_elements() != 3)
	{
		cerr << "Wrong number of annihilations after second test step." << endl;
		all_ok = false;
	}
	
	if (not (repli.rf_left.get_n_elements() == 2))
	{
		cerr << "Wrong number of left to right moving forks after second test step." << endl;
		all_ok = false;
	}
	
	if (not (repli.rf_right.get_n_elements() == 2))
	{
		cerr << "Wrong number of right to left moving forks after second test step." << endl;
		all_ok = false;
	}
	
	return all_ok;
}

int main(int argc, char* argv[])
{
	return step_test();
}
