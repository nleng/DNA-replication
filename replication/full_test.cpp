
#include <iostream>
#include <cmath>
#include "replicator.hh"

using std::cerr;
using std::endl;
using std::fabs;

bool full_test()
{
	bool all_ok = true;
	
	parameter_set params;
	params.strand_length = 1.e7;
	params.n_origins = 1e4;
	params.debug = true;
	
	chromatin_type t1;
	t1.fork_speed = 28.;
	t1.base = 0.1;
	t1.sigma = 0.5e6;
	chromatin_type t2;
	t2.fork_speed = 28.;
	t2.base = 0.01;
	t2.sigma = 0.5e6;
	params.chromatins.push_back(t1);
	params.chromatins.push_back(t2);
	
	double patch_size = 5.e6;
	size_t n_chromosomes = 10;
	double chromlim = 0.;
	double curbound = 0.;
	size_t cur_type = 0.;
	while (curbound < params.strand_length)
	{
		boundary b1(curbound, false, cur_type, cur_type==0?1:0);
		if (curbound >= chromlim)
		{
			b1.is_chromosome_boundary = true;
			chromlim += params.strand_length/n_chromosomes;
		}
		params.boundaries.push_back(b1);
		if (cur_type)
		{
			cur_type = 0;
		}
		else
		{
			cur_type = 1;
		}
		curbound += patch_size;
	}
	
	params.capacity = 5000;
	params.rate = 0.00092143193;
	
	params.write_to_disk = true;
	replicator repli;
	
	repli.set_parameters(params);
	
	repli.run();
	
	cout << "Jepp" << endl;
	/*
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
	
	if (not (fabs(a2.endtime - (0.5*repli.strand_length/t2.fork_speed + repli.current_time)) < 1e-10))
	{
		cerr << "Wrong time for annihilation of replaced fork." << endl;
		all_ok = false;
	}
	
	origin o2(0.2*repli.strand_length, 0);
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
	*/
	return all_ok;
}

int main(int argc, char* argv[])
{
	return full_test();
}
