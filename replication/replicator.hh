
#ifndef REPLICATOR_HH
#define REPLICATOR_HH


#include "simtools/binary_heap.hh"
#include "simtools/red_black_tree.hh"
#include "simtools/exceptions.hh"
#include "fork.hh"
#include "annihilation.hh"
#include "boundary.hh"
#include "origin.hh"
#include "parameter_set.hh"
#include <utility>
#include "chromatin_type.hh"
#include <vector>
#include <set>
#include <fstream>
#include <tr1/random>
#include <boost/generator_iterator.hpp>

using std::vector;
using std::set;
using std::pair;
using std::ofstream;
using simtools::binary_heap;
using simtools::red_black_tree;
using simtools::RuntimeError;

typedef std::tr1::mt19937 engine_t;
typedef std::tr1::uniform_int<size_t> uniform_int_t;
typedef std::tr1::uniform_real<long double> uniform_real_t;
typedef std::tr1::variate_generator<std::tr1::mt19937&, uniform_int_t> uni_int_gen;
typedef std::tr1::variate_generator<std::tr1::mt19937&, std::tr1::uniform_real<long double> > vg_real_t;
typedef boost::generator_iterator<std::tr1::variate_generator<std::tr1::mt19937&, std::tr1::uniform_real<long double> > > real_it_t;

bool check_anni(const binary_heap<annihilation> &inanni);
bool check_anni_pointers(const binary_heap<annihilation> &inanni);

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
	if (not all_ok)
		throw RuntimeError("Red black tree inconsistency detected.");
	return all_ok;
}

class replicator
{
	public:
		replicator();
		replicator(const replicator &right);
		
		void run();
		void set_parameters(const parameter_set &params);
		
		friend bool test_replicator();
		friend bool firing_test();
		friend bool step_test();
		
		void clear();
		
		bool test_annihilation_consistency();
// 		friend bool full_test();
	protected:
		red_black_tree<origin> origins;
		red_black_tree<fork_t> rf_left;		// Forks moving from left to right
		red_black_tree<fork_t> rf_right;	// Forks moving from right to left
		
		red_black_tree<boundary> boundaries;
		
		void clean_origins(double left, double right);
		pair<double, int> origin_prob(red_black_tree<origin>::const_iterator p) const;
		bool check_skip(red_black_tree<origin>::const_iterator p) const;
		void fire_origin(red_black_tree<origin>::iterator tofire, unsigned int inid);
		pair<red_black_tree<origin>::iterator, int> get_firing_origin();
		void fill_by_firing();
		
		void calc_firing_by_limiter_growth();
		void fill_origins();
		void checkall_annihilations() const;
		void store_lifetime(const double lifetime, const size_t chromatin);
		
		void push_left_anni(fork_t &infork, bool preserve_id);
		void push_right_anni(fork_t &infork, bool preserve_id);
		
		void step();
		void pre_output();
		void output(double out_time);
		void post_output();
		pair<vector<long double>, vector<unsigned int> > get_covered_origins() const;
		void fork_consistency_check() const;
		unsigned int count_origins(double start, double end);
		unsigned int count_forks(double start, double end);
		
		unsigned int fork_counter;
		unsigned int anni_counter;
		unsigned int cascade_counter;
		unsigned int step_counter;
		long double current_time;
		long double last_output;
		
		parameter_set parameters;
		
		binary_heap<double> new_pairs;
		binary_heap<annihilation> annihilations;
		binary_heap<change_parameter> parameter_changes;
		vector<chromatin_type> chromatins;
		// nicor
		vector<bool> chromosomeStarted;
		double max_induced_firing;
		// rocin
		vector<unsigned int> trig_firing; // The number of origins that fire triggered in each chromatin type.
		vector<unsigned int> untrig_firing; // Same for untriggered ones.
		vector<unsigned int> anni_count; // Total count of annihilations in each chromatin type.
		vector< vector<unsigned int> > *lifetime_count;	// Counters for fork lifetime bins. Each entry is a vector containing one value for each chromatin type.
		list<pair<double, double> > covered_at_time; // position/time pairs for the note_when_covered points.
		
		engine_t rng_engine;
		uniform_real_t rng_real;
		vg_real_t VG_real;
		real_it_t real_it;
		
		long double total_replicated;
//		double gaussProb = 0.05;
		unsigned int n_fired_origins;;
		
		size_t current_limit; // Current value of the limiter.
		size_t n_active_pairs; // Number of curretly active fork pairs.
		
		vector<unsigned int> n_crossovers; // How many crossovers happened between any two chromatin types.
		
		set<double> firing_positions;
		
		ofstream rate_file;	// Replication rate and more go into this file.
		ofstream co_file; 	// Covered origins file.
		ofstream anni_file;
		ofstream ip_file;
		ofstream lifetime_file;
		ofstream ori_d_file;
		ofstream fdd_file;
		ofstream fork_file;
		// nicor binary
		FILE* fork_file_binary;
		// rocin
		ofstream log_file;
		
  		unsigned int runcounter;
		bool first_fired;
		bool simulation_done;
		bool firing_prob_summation_mode;
};

#endif // REPLICATOR_HH
