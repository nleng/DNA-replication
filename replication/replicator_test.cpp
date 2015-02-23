
#include <iostream>
#include <cmath>
#include "replicator.hh"

using std::cerr;
using std::endl;
using std::fabs;

bool test_replicator()
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
	
	repli.fire_origin(p1, 0);
	
	if (repli.annihilations.get_n_elements() != 2)
	{
		cerr << "Wrong number of forks after first firing event." << endl;
		all_ok = false;
	}
	
	annihilation a1 = repli.annihilations.get_value(0);
	annihilation a2 = repli.annihilations.get_value(1);
	if (not((a1.at == SOLO_REMOVAL) and (a1.left)))
	{
		cerr << "Wrong annihilation type for second fork." << endl;
		all_ok = false;
	}
	
	if (not (fabs(a1.endtime - 0.8*repli.parameters.strand_length/t1.fork_speed) < 1e-10))
	{
		cerr << "Wrong end time for second fork annihilation." << endl;
		all_ok = false;
	}
	
	if (not((a2.at == SOLO_REMOVAL) and (not a2.left)))
	{
		cerr << "Wrong annihilation type for first fork." << endl;
		all_ok = false;
	}
	
	if (not (fabs(a2.endtime - 0.2*repli.parameters.strand_length/t1.fork_speed) < 1e-10))
	{
		cerr << "Wrong end time for first fork annihilation." << endl;
		all_ok = false;
	}
	
	repli.current_time = 2.;
	
	origin o2(0.3*repli.parameters.strand_length, 0);
	red_black_tree<origin>::iterator p2 = repli.origins.insert(o2);
	
	repli.fire_origin(p2, 0);
	
	if ((repli.annihilations.get_n_elements() != 3) or (repli.rf_left.get_n_elements() != 2) or (repli.rf_right.get_n_elements() != 2))
	{
		cerr << "Wrong number of forks after second firing event." << endl;
		all_ok = false;
	}
	
	annihilation a4 = repli.annihilations.get_value(1);
	annihilation a5 = repli.annihilations.get_value(2);
	annihilation a6 = repli.annihilations.get_value(3);
	
	if (not((a4.at == SOLO_REMOVAL) and (not a4.left)))
	{
		cerr << "Wrong annihilation type for second fork." << endl;
		all_ok = false;
	}
	
	if (a6.at != FORK_FORK)
	{
		cerr << "Wrong annihilation type for first and fourth fork." << endl;
		all_ok = false;
	}
	
	if (not((a5.at == SOLO_REMOVAL) and (a5.left)))
	{
		cerr << "Wrong annihilation type for third fork." << endl;
		all_ok = false;
	}
	
	if (not (fabs(a4.endtime - 0.2*repli.parameters.strand_length/t1.fork_speed) < 1e-10))
	{
		cerr << "Wrong end time for second fork annihilation." << endl;
		all_ok = false;
	}
	
	if (not (fabs(a5.endtime - (2. +0.7*repli.parameters.strand_length/t1.fork_speed)) < 1e-10))
	{
		cerr << "Wrong end time for second fork annihilation." << endl;
		all_ok = false;
	}
	
	if (not (fabs(a6.endtime - (2.+ (0.1*repli.parameters.strand_length -2.*t1.fork_speed)/2./t1.fork_speed)) < 1e-10))
	{
		cerr << "Wrong end time for first and fourth fork annihilation." << endl;
		all_ok = false;
	}
	
	red_black_tree<fork_t>::iterator rbtp = repli.rf_left.find_biggest(0);
	if (rbtp->value.anni != 2)
	{
		cerr << "Third fork has annihilation inconsistency." << endl;
		all_ok = false;
	}
	
	rbtp = repli.rf_left.find_biggest(1);
	if (rbtp->value.anni != 3)
	{
		cerr << "First fork has annihilation inconsistency." << endl;
		all_ok = false;
	}
	
	rbtp = repli.rf_right.find_biggest(0);
	if (rbtp->value.anni != 3)
	{
		cerr << "Fourth fork has annihilation inconsistency." << endl;
		all_ok = false;
	}
	
	rbtp = repli.rf_right.find_biggest(1);
	if (rbtp->value.anni != 1)
	{
		cerr << "Second fork has annihilation inconsistency." << endl;
		all_ok = false;
	}
	
	return all_ok;
}

int main(int argc, char* argv[])
{
	return test_replicator();
}
