
#include "parameter_set.hh"

parameter_set::parameter_set()
	: strand_length(1.e7), n_origins(1e4), capacity(5000), rate(40.), inhibition_distance(1e3), output_delta(100), lifetime_bin_width(100.), firing_distance_bin_width(100.), rate_file_name("res_rate_data.txt"), covered_origin_file_name("res_covered_origins.txt"), ori_file_name("res_origins.txt"), anni_file_name("res_annihilations.txt"), ip_file_name("res_induced_patches.txt"), lifetimes_file_name("res_lifetimes.txt"), ori_d_file_name("res_debug_origins.txt"), fdd_file_name("res_fdd.txt"), fork_file_name("res_forks.txt"), log_file_name("replication_log.txt"), note_when_covered_file_name("red_nwc.txt"), use_inhibition(false), use_induced_limit(false), use_origin_list(false), use_alternate_induced_firing(false), write_to_disk(false), write_covered_origins(false), write_induced_patches(false), write_fdd(false), write_note_when_covered(false), dump_fork_positions(false), insane_origin_debug_level(false), debug(false)
{
}
	
parameter_set& parameter_set::operator= (const parameter_set &right)
{
	// nicor
	chromosomeStarted = right.chromosomeStarted;
	max_induced_firing = right.max_induced_firing;
	// rocin
	strand_length = right.strand_length;
	n_origins = right.n_origins;
	boundaries = right.boundaries;
	chromatins = right.chromatins;
	origins = right.origins;
	change_params = right.change_params;
	earlyfire = right.earlyfire;
	note_when_covered = right.note_when_covered;
	capacity = right.capacity;
	rate = right.rate;
	inhibition_distance = right.inhibition_distance;
	output_delta = right.output_delta;
	lifetime_bin_width = right.lifetime_bin_width;
	firing_distance_bin_width = right.firing_distance_bin_width;
	rate_file_name = right.rate_file_name;
	covered_origin_file_name = right.covered_origin_file_name;
	ori_file_name = right.ori_file_name;
	anni_file_name = right.anni_file_name;
	ip_file_name = right.ip_file_name;
	lifetimes_file_name = right.lifetimes_file_name;
	ori_d_file_name = right.ori_d_file_name;
	fdd_file_name = right.fdd_file_name;
	fork_file_name = right.fork_file_name;
	log_file_name = right.log_file_name;
	note_when_covered_file_name = right.note_when_covered_file_name;
	use_inhibition = right.use_inhibition;
	use_induced_limit = right.use_induced_limit;
	use_origin_list = right.use_origin_list;
	use_alternate_induced_firing = right.use_alternate_induced_firing;
	write_to_disk = right.write_to_disk;
	write_covered_origins = right.write_covered_origins;
	write_induced_patches = right.write_induced_patches;
	write_fdd = right.write_fdd;
	write_note_when_covered = right.write_note_when_covered;
	dump_fork_positions = right.dump_fork_positions;
	insane_origin_debug_level = right.insane_origin_debug_level;
	debug = right.debug;
	
	return *this;
}

#ifdef USES_BOOST_PYTHON
void  paremeter_set::set_boundaries(const boost::python::list &inlist)
{
	boundaries.clear();
		// Takes a list of lists with the inner lists containing a double (position) and size_t (index) each.
	size_t n_bs = boost::python::len(inlist);
		
	for(size_t i = 0; i < n_bs; i++)
	{
		boost::python::list templist = boost::python::extract<boost::python::list>(inlist[i]);
		if (boost::python::len(templist) < 4)
		{
			throw RuntimeError("List entry too short in boundary list. Must have at least 4 entries.");
		}
		else
		{
			boundaries.insert( boundary(
					   boost::python::extract<double>(templist[0]), 
					boost::python::extract<bool>(templist[1]),
					boost::python::extract<size_t>(templist[2]),
					boost::python::extract<size_t>(templist[3])));
		}	
	}
}
	
boost::python::list parameter_set::get_boundaries() const
{
	boost::python::list outlist;
	for(set< boundary >::iterator p = boundaries.begin(); p != boundaries.end(); ++p)
	{
		boost::python::list tmplist;
		tmplist.append(boost::python::object(p->position));
		tmplist.append(boost::python::object(p->is_chromosome_boundary));
		tmplist.append(boost::python::object(p->left));
		tmplist.append(boost::python::object(p->right));
		outlist.append(tmplist);
	}
	return outlist;
}
	
void parameter_set::set_chromatins(const boost::python::list &inchroms)
{
	chromatins.clear();
		// Takes a list of chromatin_zones.
	size_t n_bs = boost::python::len(inchroms);
	for(size_t i = 0; i < n_bs; i++)
	{
		chromatins.push_back(boost::python::extract<chromatin_type>(inchroms[i]));
	}
}
	
boost::python::list parameter_set::get_chromatins() const
{
	boost::python::list outlist;
	for(vector<chromatin_type>::const_iterator p = chromatins.begin(); p != chromatins.end(); ++p)
	{
		outlist.append(*p);
	}
	return outlist;
}
#endif
