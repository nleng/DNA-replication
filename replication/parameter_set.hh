
#ifndef PARAMETER_SET
#define PARAMETER_SET

#ifdef USES_BOOST_PYTHON
#include <boost/python.hpp>
#endif

#include <string>
#include <vector>
#include <set>
#include "boundary.hh"
#include "chromatin_type.hh"
#include "change_parameter.hh"

using std::string;
using std::vector;
using std::set;

struct parameter_set
{
	parameter_set();
	
	parameter_set& operator= (const parameter_set &right);
	
	// Chromatin data:
	double strand_length;
	unsigned int n_origins;
	
	vector<boundary> boundaries;
	vector<chromatin_type> chromatins;
	vector<double> origins;
	// nicor
	vector<bool> chromosomeStarted;
	double max_induced_firing;
	// rocin
	vector<change_parameter> change_params;
	set<unsigned> earlyfire; // Indices of origins that fire early (i.e. at time 0)
	set<double> note_when_covered; // Positions of which the replication time should be noted.
	
	// Limiter parameters:
	double capacity;
	double rate;
	
	// If neccessary: Firing inhibition parameters:
	double inhibition_distance;
	
	// Technical simulation parameters:
	double output_delta;
	double lifetime_bin_width;
	
	double firing_distance_bin_width;
	
	string rate_file_name;
	string covered_origin_file_name;
	string ori_file_name;
	string anni_file_name;
	string ip_file_name;
	string lifetimes_file_name;
	string ori_d_file_name;
	string fdd_file_name;	// Firing distance distribution.
	string fork_file_name;
	string log_file_name;
	string note_when_covered_file_name;
	
	bool use_inhibition;
	bool use_induced_limit; // This means that induced firing only happens when the induced probability is bigger than a lower limit.
	bool use_origin_list; // If this parameter is set, don't generate randomly distributed origins but choose the contents of "origins" instead.
	bool use_alternate_induced_firing;
	bool write_to_disk;
	bool write_covered_origins;
	bool write_induced_patches; // If true, the number of contiguous patches of induced firing will be written to disk.
	bool write_fdd; // Firing distance distribution.
	bool write_note_when_covered;
	bool dump_fork_positions;	// If set, write all fork positions to disk.
	bool insane_origin_debug_level;
        bool debug;
	
#ifdef USES_BOOST_PYTHON
	void set_boundaries(const boost::python::list &inlist);
	boost::python::list get_boundaries() const;
	void set_chromatins(const boost::python::list &inchroms);
	boost::python::list get_chromatins() const;
#endif // USES_BOOST_PYTHON
	
};

#endif //PARAMETER_SET
