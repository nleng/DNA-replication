#include <iostream>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <set>
#include <list>
#include <utility>
#include "replicator.hh"
#include "simtools/exceptions.hh"
#include <boost/math/special_functions/gamma.hpp>

using std::cerr;
using std::endl;
using std::set;
using std::vector;
using std::fabs;
using std::string;
using simtools::tostr;
using std::pow;
using std::exp;
using boost::math::tgamma;
using std::ifstream;
using std::istringstream;
using std::list;
using std::pair;
using std::scientific;

typedef std::tr1::normal_distribution<long double> normal_t;
typedef std::tr1::variate_generator<std::tr1::mt19937&, std::tr1::normal_distribution<long double> > vg_norm_t;
typedef boost::generator_iterator<std::tr1::variate_generator<std::tr1::mt19937&, std::tr1::normal_distribution<long double> > > real_norm_it_t;

bool run()
{
	bool all_ok = true;
	
	unsigned n_origins = 500000;
	double inhibition_distance = 55000.;
	double sigma = 2.4e5;			// Sigma of induced firing gaussian
	double induced_limit = 0.0;		// Cutoff probability density for induced firing. 
	size_t n_runs = 1;
	double limiter = 6000;			// This is the number of fork _pairs_.
	double initial_slope = 0.001111111111; 	// This parameter determines the initial limiter growth.	entspricht 1/(15 min), aus der funktion Lmax*(1-exp(-t/tau))
	double fork_speed = 28.;

	// Load the chromosome data that was taken from the USC genome browser. UCSC
	ifstream infile("chromosome_patchlists_inactiveX.txt");
       	vector<pair<size_t, list<double> > > band_data;
	
	while(infile)
	{
		list<double> t_bands;
		size_t chromtype = 0;
		string t_line;
		getline(infile, t_line);
		istringstream line(t_line);
		
		line.ignore(2,'\n');
		line >> chromtype;
		
		while(not line.eof())
		{
			double tmp;
			line >> tmp;
			t_bands.push_back(tmp);
		}
		
		pair<size_t, list<double> > rpair(chromtype, t_bands);
		band_data.push_back(rpair);
	}

	ifstream infile2("chromatinType.txt");
    vector<vector<int>> chromType_data;
	while(infile2)
	{
		vector<int> t_bands;
		string t_line;
		getline(infile2, t_line);
		istringstream line(t_line);

		while(not line.eof())
		{
			int tmp;
			line >> tmp;
			t_bands.push_back(tmp);
		}

		chromType_data.push_back(t_bands);
	}


	// hier inactiveX chromosome
	vector<size_t> n_chromosomes = {5,3,3,2,6,3,5,3,5,3,3,3,3,3,3,3,4,2,3,3,3,3,1,1};
	
	set<double> tnwc;
	ifstream nwcinfile("nwc_positions.txt");
	while(nwcinfile)
	{
		string t_line;
		getline(nwcinfile, t_line);
		istringstream line(t_line);
		
		double tpos = 0.;
		line >> tpos;
		if (line != 0)
			tnwc.insert(tpos);
	}
	
	// Output directory.
	string dirname("/home/corni/aaReplication/replicationSim/resultsEclipse/paper_hela_realistic_forks/");
	parameter_set params;
	
	params.write_covered_origins = false;
	params.debug = false;
	params.dump_fork_positions = true;
	params.use_alternate_induced_firing = false;
	params.write_note_when_covered = true;
	params.insane_origin_debug_level = false;
	
	replicator repli;
	
 	params.capacity = limiter;
	params.rate = initial_slope;
	params.write_to_disk = true;
	params.write_fdd = true;
	params.write_induced_patches = true;
	params.firing_distance_bin_width = 10000;
	
	engine_t rng_engine(static_cast<std::tr1::mt19937::result_type>(42));
	
	normal_t rngr(0., 1.);
	vg_norm_t VGr(rng_engine, rngr);
	real_norm_it_t rit(&VGr);
	
	uniform_real_t rng_real(0.,1.);
	vg_real_t VG_real(rng_engine, rng_real);
	real_it_t real_it(&VG_real);
	
	for (size_t i = 0; i < n_runs; i++)
	{
		params.boundaries.clear();
		params.chromatins.clear();
		params.origins.clear();
		params.change_params.clear();
		params.earlyfire.clear();
			
		params.n_origins = n_origins;
				
		chromatin_type t1;
		t1.fork_speed = 6.;
		t1.base = 0.8;	// spontaneous firing probability for euchromatin
		t1.sigma = sigma;
		t1.induced_limit = induced_limit;
		chromatin_type t2;
		t2.fork_speed = 6.;
		t2.base = 0.0;	// spontaneous firing probability for heterochromatin
      	t2.sigma = sigma;
		t2.induced_limit = induced_limit;
		chromatin_type t3;
		t3.fork_speed = 6.;
		t3.base = 0.05;	// spontaneous firing probability for facultative heterochromatin
      	t3.sigma = sigma;
		t3.induced_limit = induced_limit;
		
		params.chromatins.clear();
		params.chromatins.push_back(t1);
		params.chromatins.push_back(t2);
		params.chromatins.push_back(t3);
		
		// hier ?????
		params.max_induced_firing = 1.0;
		params.chromosomeStarted = {};
//		for(int a=0;a<n_chromosomes.size();a++){
//			for(int b=0;b<n_chromosomes[a];b++){
//				params.chromosomeStarted.push_back(false);
//			}
//		}
		// ?????

		// alle speed aenderungen werden eingetragen
		change_parameter cp1;
		cp1.command = CHANGE_REPLI_SPEED;

		for (int vv = 1; vv <= 14; vv++)
		{
			double time = (double) vv*720.;
			cp1.time = time;
			cp1.value[0] = 6.+(fork_speed-6.)*time/10080.;
			cp1.chromatin = 0;
			params.change_params.push_back(cp1);
			cp1.value[0] = 6.+(fork_speed-6.)*time/10080.;
			cp1.chromatin = 1;
			params.change_params.push_back(cp1);
			cp1.value[0] = 6.+(fork_speed-6.)*time/10080.; // parameter[9]/4.*(1.+3.*time/10800.);
			cp1.chromatin = 2;
			params.change_params.push_back(cp1);
		}
		
		params.use_induced_limit = true;
		params.use_inhibition = true;
		params.inhibition_distance = inhibition_distance;
		
		params.use_origin_list = true;
		params.origins.reserve(n_origins);
		
		params.boundaries.clear();
				
		double eu_mass = 0.;
		double het_mass = 0.;
		double fac_mass = 0;
		
		double curbound = 0.;

//		size_t cur_type = band_data[0].first;
		size_t cur_type = chromType_data[0][0];
		size_t total_chromo_counter = 0;


		boundary b1(curbound, true, chromType_data[0][1], cur_type, -1, total_chromo_counter);
		params.boundaries.push_back(b1);
		for (size_t jj = 0; jj < n_chromosomes.size(); ++jj)
		{
			for (size_t ii = 0; ii < n_chromosomes[jj]; ++ii)
			{
				size_t bandCounter = 1;
//				cur_type = band_data[jj].first;	// daniel
				cur_type = chromType_data[jj][0];	// nicor
//				cur_type = chromType_data[total_chromo_counter][0];	// nicor
				params.boundaries[params.boundaries.size()-1].right = cur_type;

				for (list<double>::iterator q = band_data[jj].second.begin(); q != band_data[jj].second.end(); ++q)
				{
					// das hier sind die chomatin type boundaries
					size_t last_type = cur_type;
//					cur_type = cur_type==0?1:0;	// daniel
					// jj durch total_chromo_counter ersetzen, falls facultative aus 3D position, insgesamt an 3 stellen
					bandCounter = bandCounter%chromType_data[jj].size();	// nicor
					cur_type = chromType_data[jj][bandCounter];	// nicor
					cout << jj << " " << bandCounter << " " << chromType_data[jj].size() << " " << last_type << " " << cur_type << endl;

					curbound += *q;
					b1.is_chromosome_boundary = false;	// hier ?????
					b1.left = last_type;
					b1.right = cur_type;
					b1.left_chromosome = total_chromo_counter;
					b1.right_chromosome = total_chromo_counter;
					b1.position = curbound;
					
					params.boundaries.push_back(b1);
					++bandCounter;
					// falls zwischen jeder zone eine chromosome_boudary ist
//					params.chromosomeStarted.push_back(false);
//					++total_chromo_counter;
				}
				// chromosome boundaries
				params.boundaries[params.boundaries.size()-1].is_chromosome_boundary = true;
				params.boundaries[params.boundaries.size()-1].right_chromosome = total_chromo_counter+1;
				// falls nur zwischen chromosomen eine chromosome_boundary ist
				params.chromosomeStarted.push_back(false);
				++total_chromo_counter;

			}
		}
		
		params.strand_length = curbound;
		params.note_when_covered = tnwc;
		
		set<pair<double, int> > temporis;
		set<double> temporis2;
		
		for (unsigned k = 0; k < n_origins; ++k)
		{
			double tmpval = params.strand_length*(*real_it++);
			
			while (temporis2.find(tmpval) != temporis2.end())
				tmpval = params.strand_length*(*real_it++);
			
			params.origins.push_back(tmpval);
			temporis.insert(pair<double,int>(tmpval,k));
			temporis2.insert(tmpval);
		}
		
		vector<boundary>::iterator q = params.boundaries.begin();
		if (params.boundaries.size() > 1)
			q++;
		else
		{
			cerr << "There are no two boundaries present. Aborting" << endl;
			break;
		}
		
		unsigned int n_eus = 0;
		unsigned int n_hets = 0;
		unsigned int n_facs = 0;
		double lastpos = 0.;
		set<boundary> bset;
		for (vector<boundary>::iterator p = params.boundaries.begin(); p != params.boundaries.end(); ++p)
		{
		    	bset.insert(*p);
		}

		for (vector<boundary>::iterator p = params.boundaries.begin(); q != params.boundaries.end(); ++p, ++q)
		{
			if ((p->right == 0) and (q->left == 0))
			{
				eu_mass += q->position - p->position;
				n_eus++;
			}
			else if ((p->right == 1) and (q->left == 1))
			{
				het_mass += q->position - p->position;
				n_hets++;
			}
			else if ((p->right == 2) and (q->left == 2))
			{
				fac_mass += q->position - p->position;
				n_facs++;
			}
			
		}
		cout << "SimNr: " << i << endl;
		cout << std::scientific << "Eu: " << eu_mass << "  Het: " << het_mass << "  Fac: " << fac_mass << endl;
		cout << std::scientific << "Avg eu size: " << eu_mass/n_eus << "  Avg het size: " << het_mass/n_hets << "  Avg fac size: " << fac_mass/n_facs << endl;
		cout << std::scientific << "N eu: " << n_eus << "  N het: " << n_hets << "  N fac: " << n_facs << endl;
		
		repli.clear();
		string identifier = tostr(limiter) + string("_") + tostr(sigma) + string("_") + tostr(inhibition_distance) + string("_") + tostr(fork_speed) + string("_") + tostr(n_origins) + string("_") + tostr(i) + string("_run");
		 // replicated amount per time
		params.rate_file_name = dirname + string("res_rates_") + identifier + string(".txt");
		// covered origins
		params.covered_origin_file_name = dirname + string("res_cov_origins_") + identifier + string(".txt");
		params.ori_file_name = dirname + string("res_origins_") + identifier + string(".txt");
		params.anni_file_name = dirname + string("res_annihilations_") + identifier + string(".txt");
		// cluster lifetimes
		params.lifetimes_file_name = dirname + string("res_lifetimes_") + identifier + string(".txt");
		params.fdd_file_name = dirname + string("res_fdd_") + identifier + string(".txt");
		params.ip_file_name = dirname + string("res_ip_") + identifier + string(".txt");
		params.fork_file_name = dirname + string("res_fork_") + identifier + string(".txt");
		// note when covered (timestap when origin is passed)
		params.note_when_covered_file_name = dirname + string("res_nwc_") + identifier + string(".txt");
		
		params.ori_d_file_name = dirname + string("res_insane_") + identifier + string(".txt");
		
		repli.set_parameters(params);
		repli.run();
	}
	
	cout << "Completed." << endl;
	return all_ok;
}

int main(int argc, char* argv[])
{
	return run();
}
