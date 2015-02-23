

#include "simtools/binary_sum_search.hh"
#include "simtools/twoe.hh"
#include "replicator.hh"
#include "exceptions.hh"
#include <cmath>
#include <map>
#include <algorithm>
#include <iomanip>

using std::exp;
using std::log;
using std::fabs;
using std::map;
using std::min;
using std::max;
using simtools::binary_sum_search;
using simtools::tostr;

replicator::replicator()
  : fork_counter(0), anni_counter(0), cascade_counter(0), step_counter(0), current_time(0.), last_output(0.), lifetime_count(NULL), rng_engine(static_cast<std::tr1::mt19937::result_type>(23)), rng_real(0., 1.), VG_real(rng_engine, rng_real), real_it(&VG_real), total_replicated(0.), n_fired_origins(0), current_limit(0.), n_active_pairs(0.), runcounter(0), first_fired(false), simulation_done(false), firing_prob_summation_mode(false)
{
}

// bool check_anni(const binary_heap<annihilation> &inanni)
// {
// 	for (map<long unsigned, long unsigned>::const_iterator p =  inanni.positions.begin(); p != inanni.positions.end(); ++p)
// 	{
// 		pair<annihilation, long unsigned int> temp = inanni.data[p->second];
// 		if (not (temp.first.first_fork->parent == 0))
// 			if(not ((temp.first.first_fork->parent->left == temp.first.first_fork->me) or (temp.first.first_fork->parent->right == temp.first.first_fork->me)))
// 	  		        throw RuntimeError("Annihilations show backlinking inconsistency.");
// 	}
// 	return true;
// }
// 
// bool check_anni_pointers(const binary_heap<annihilation> &inanni)
// {
// 	for (map<long unsigned, long unsigned>::const_iterator p =  inanni.positions.begin(); p != inanni.positions.end(); ++p)
// 	{
// 		map<long unsigned, long unsigned>::const_iterator q =  p;
// 		++q;
// 		while(q != inanni.positions.end())
// 		{
// 			if (inanni.data[p->second].first.hasequal_fork(inanni.data[q->second].first))
// 			        throw RuntimeError("A fork appears in more than one annihilation.");
// 			++q;
// 		}
// 	}
// 	return true;
// }

void replicator::clear()
{
	origins.clear();
	rf_left.clear();
	rf_right.clear();
	boundaries.clear();
	
	fork_counter = 0;
	anni_counter = 0;
	cascade_counter = 0;
	step_counter = 0;
	new_pairs.clear();
	annihilations.clear();
	parameter_changes.clear();
	chromatins.clear();
	trig_firing.clear();
	untrig_firing.clear();
	anni_count.clear();
	if (lifetime_count)
	{
		lifetime_count->clear();
		delete lifetime_count;
	}
	lifetime_count = NULL;
	
	current_time = 0.;
	last_output = 0.;
	total_replicated = 0.;
	n_fired_origins = 0;
	current_limit = 0;
	n_active_pairs = 0;
	
	n_crossovers.clear();
	firing_positions.clear();
	covered_at_time.clear();
	
	first_fired = false;
	simulation_done = false;
	firing_prob_summation_mode = false;
}

void replicator::run()
{
	calc_firing_by_limiter_growth();
	if (boundaries.get_n_elements() < 2)
	{
		throw RuntimeError("Simulation cannot start without at least two boundaries.");
	}
	
	pre_output();
	output(current_time);
	 
	while(((annihilations.get_n_elements() > 0) or (new_pairs.get_n_elements() > 0)) and (not simulation_done))
	{
		step();
	}
	post_output();
	runcounter++;
}

void replicator::step()
{
        if (parameters.debug)
	{
	        checkall_annihilations();
//                 check_anni(annihilations);
	}
	
	// Perform the next step of this simulation.
	
	// Determine, which event is the next one to trigger action.
	long double nexttime = 1.e20;
	bool isanni = false;
	bool isparamchange = false;
	
	// abfrage ob neue paare zur verfuegung stehen
	if (new_pairs.get_n_elements() > 0)
	{
		pair<double, unsigned long> resp = new_pairs.get_smallest();
		nexttime = resp.first;
	}
	
	// abfrage ob annihilation ansteht
	if (annihilations.get_n_elements() > 0)
	{
		pair<annihilation, unsigned long> resp = annihilations.get_smallest();
		if (resp.first.endtime < nexttime)
		{
			nexttime = resp.first.endtime;
			isanni = true;
		}
	}
	
	// abfrage ob parameter geaendert wurden (geschwindigkeit)?
	if (parameter_changes.get_n_elements() > 0)
	{
		pair<change_parameter, unsigned long> resp = parameter_changes.get_smallest();
		if (resp.first.time < nexttime)
		{
		  	isanni = false;
			isparamchange = true;
			nexttime = resp.first.time;
		}
	}
	
	if (nexttime < current_time)
	{
		throw RuntimeError("Time ordering inconsistency.");
	}
	
	// aktualisiert die zeit? output() funktion noch einmal anschauen ????? 
	while(last_output + parameters.output_delta < nexttime)
	{
		last_output += parameters.output_delta;
		output(last_output);
	}
	
	// eigentlicher annihilations-event, der binary_heap mit allen funktionen, wie .remove_smallest() ist in den simtools definiert
	if (isanni)
	{
		pair<annihilation, unsigned long> resp = annihilations.remove_smallest();
		
		// bestimmung von chromatin-type
		size_t chrom_type = resp.first.first_fork->chromatin;

		// FORK_FORK ist der annihilations-typ, bei den zwei forks zusammenstossen
		if (resp.first.at == FORK_FORK)
		{
			anni_count[chrom_type]++;
			anni_count[chrom_type]++;
		}
		else if (resp.first.at == SOLO_REMOVAL)
		{
			anni_count[chrom_type]++;
		}
		// fork+fork annihilation, .at gibt annihilation-typ an
		current_time = resp.first.endtime;
//		std::cout << current_time << endl;
		if (resp.first.at == FORK_FORK)
		{
			store_lifetime(current_time-resp.first.first_fork->start_time, resp.first.first_fork->chromatin);
			store_lifetime(current_time-resp.first.second_fork->start_time, resp.first.second_fork->chromatin);
			// testet, ob ?????
			if (resp.first.left)
			{
				total_replicated += resp.first.second_fork->position -resp.first.first_fork->position;
				clean_origins(resp.first.first_fork->position, resp.first.second_fork->position);
				rf_left.remove(resp.first.first_fork);
				rf_right.remove(resp.first.second_fork);
			}
			else
			{
				total_replicated += resp.first.first_fork->position - resp.first.second_fork->position;
				clean_origins(resp.first.second_fork->position, resp.first.first_fork->position);
				rf_right.remove(resp.first.first_fork);
				rf_left.remove(resp.first.second_fork);
			}
			
			n_active_pairs--;
			fill_by_firing();
		}
		// falls eine fork ueber den rand hinauslauft, hier nur chromatin zone boundaries, transition, crossover
		else if (resp.first.at == FORK_BOUNDARY)
		{
			fork_t newin(resp.first.other_boundary->position, current_time, resp.first.first_fork->origin_id, resp.first.first_fork->id, 0, 0);
			newin.start_time = resp.first.first_fork->start_time;
			
			if (resp.first.left)
			{
				size_t ncind =  resp.first.first_fork->chromatin*chromatins.size() + resp.first.other_boundary->right;
				n_crossovers[ncind]++;
				
				total_replicated += resp.first.other_boundary->position - resp.first.first_fork->position;
				clean_origins(resp.first.first_fork->position, resp.first.other_boundary->position);
				rf_left.remove(resp.first.first_fork);
				newin.chromatin = resp.first.other_boundary->right;
				newin.chromosome = resp.first.other_boundary->right_chromosome;
				push_left_anni(newin, true);
			}
			else
			{
				size_t ncind =  resp.first.first_fork->chromatin*chromatins.size() + resp.first.other_boundary->left;
				n_crossovers[ncind]++;
				
				total_replicated += resp.first.first_fork->position - resp.first.other_boundary->position;
				clean_origins(resp.first.other_boundary->position, resp.first.first_fork->position);
				rf_right.remove(resp.first.first_fork);
				newin.chromatin = resp.first.other_boundary->left;
				newin.chromosome = resp.first.other_boundary->left_chromosome;
				push_right_anni(newin, true);
			}
		}
		// falls ein chromosom verlassen wird
		else if (resp.first.at == SOLO_REMOVAL)
		{
			store_lifetime(current_time-resp.first.first_fork->start_time, resp.first.first_fork->chromatin);
			if (resp.first.left)
			{
				total_replicated += resp.first.other_boundary->position - resp.first.first_fork->position;
				clean_origins(resp.first.first_fork->position, resp.first.other_boundary->position);
				rf_left.remove(resp.first.first_fork);
			}
			else
			{
				total_replicated += resp.first.first_fork->position - resp.first.other_boundary->position;
				clean_origins(resp.first.other_boundary->position, resp.first.first_fork->position);
				rf_right.remove(resp.first.first_fork);
			}
			
			if (not ((rf_left.get_n_elements() + rf_right.get_n_elements())%2))
			{
				n_active_pairs--;
				fill_by_firing();
			}
		}
	}
	else if (not isparamchange)
	{
		pair<double, unsigned long> resp = new_pairs.remove_smallest();
		
		current_time = resp.first;
		
		if (origins.get_n_elements() > 0)
		{
			current_limit++;
			fill_by_firing();
		}
		else
		{
			if (annihilations.get_n_elements() == 0)
			{
				// End the simulation because all that is left are limiter increase events.
				simulation_done = true;
			}
		}
	}
	// falls sich parameter geaendert haben, hier eigentlich nur die geschwindigkeit oda?
	if (isparamchange)
	{
		pair<change_parameter, unsigned long> resp = parameter_changes.remove_smallest();
		
		current_time = resp.first.time;
		if (resp.first.command == CHANGE_REPLI_SPEED)
		{
			// If there was a speedchange, all forks have to be updated by revising their starting
			// positions and starting times. If we don't do this, there will be annihilations that
			// happen in the past, making the algorithm crash.
			
			// We use the current time as the switch time, not the switch time that was determined
			// earlier. The reason for this is that we want to avoid annihilations in the past.
			// Since the inaccuracy that is incurred is only in the range of tens of seconds, we
			// are willing to accept this.
			
			for (size_t i = 0; i < rf_left.get_n_elements(); i++)
			{
				red_black_tree<fork_t>::iterator p = rf_left.find_biggest(i);

				double opos = p->position;
				p->position += chromatins[p->chromatin].fork_speed*(current_time-p->time);
				p->time = current_time;
				clean_origins(opos, p->position);
			}

			for (size_t j = 0; j < rf_right.get_n_elements(); j++)
			{
				red_black_tree<fork_t>::iterator p = rf_right.find_biggest(j);

				double opos = p->position;
				p->position -= chromatins[p->chromatin].fork_speed*(current_time-p->time);
				p->time = current_time;
				clean_origins(p->position, opos);
			}
			
			// Now update the speed.
			chromatins[resp.first.chromatin].fork_speed = resp.first.value[0];
			
			// Also update all the annihilations.
			binary_heap<annihilation> tmpani = annihilations;
			annihilations.clear();
			
			while (tmpani.get_n_elements() > 0)
			{
				pair<annihilation, long unsigned> tvl = tmpani.remove_smallest();
				
				if (tvl.first.at == FORK_FORK)
				{
					tvl.first.endtime = fabs(tvl.first.first_fork->position - tvl.first.second_fork->position)/2./chromatins[tvl.first.first_fork->chromatin].fork_speed + current_time;
					
				}
				else
				{
					tvl.first.endtime = fabs(tvl.first.first_fork->position - tvl.first.other_boundary->position)/chromatins[tvl.first.first_fork->chromatin].fork_speed + current_time;
				}
				annihilations.insert_value(tvl.first, tvl.second);
			}
			
			if (not annihilations.test_consistency())
				throw RuntimeError("Annihilations inconsistent after speed change.");
			if (not rf_left.check_sortedness().first)
				throw RuntimeError("Left forks aren't sorted after speed change");
			if (not rf_right.check_sortedness().first)
				throw RuntimeError("Right forks aren't sorted after speed change");
		}
		else if (resp.first.command == CHANGE_BASE)
		{
			current_time = resp.first.time;
			
			chromatins[resp.first.chromatin].base = resp.first.value[0];
		}
		else if (resp.first.command == CHANGE_SIGMA)
		{
			current_time = resp.first.time;
			
			chromatins[resp.first.chromatin].sigma = resp.first.value[0];
		}
	}
	
	step_counter++;
}

void replicator::fork_consistency_check() const
{
	// Test if the current number of forks is consistent with the number of annihilations, crossovers and so on.
	
	vector<unsigned int> forkcnt;
	for (size_t i = 0; i < chromatins.size(); i++)
	{
		forkcnt.push_back(0);
	}
	
	for (size_t j = 0; j < rf_left.get_n_elements(); j++)
	{
		red_black_tree<fork_t>::const_iterator p = rf_left.find_biggest_const(j);
		forkcnt[p->chromatin]++;
	}
	
	for (size_t j = 0; j < rf_right.get_n_elements(); j++)
	{
		red_black_tree<fork_t>::const_iterator p = rf_right.find_biggest_const(j);
		forkcnt[p->chromatin]++;
	}
	
	for (size_t k = 0; k < chromatins.size(); k++)
	{
		unsigned int moveout = 0;
		for (size_t ii = k*(chromatins.size()); ii < (k+1)*(chromatins.size()); ii++)
		{
			moveout += n_crossovers[ii];
		}
		
		unsigned int movein = 0;
		for (size_t jj = k; jj < chromatins.size()*chromatins.size(); jj += chromatins.size())
		{
			movein += n_crossovers[jj];
		}
		
		if (2*trig_firing[k] + 2*untrig_firing[k] + movein - moveout - anni_count[k] != forkcnt[k])
		{
			throw RuntimeError(string("Inconsistent fork number for chromatin ") + tostr(k) + string(" is ") + tostr(forkcnt[k]) + string(" but should be ") + tostr(2*trig_firing[k] + 2*untrig_firing[k] + movein - moveout - anni_count[k]) + string("."));
		}
	}
}

void replicator::clean_origins(double left, double right)
{
	// Remove all origins that exist between left and right. We have to do this to avoid firing of origins in areas that already have been replicated. (Especially complete cromosomes)
	
	if (origins.get_n_elements() == 0)
		return;
	origin otmp(left,0,-1);
	size_t cntr = 0;
	red_black_tree<origin>::iterator p = origins.find_same(otmp);
	if (p == NULL)
		p = origins.find_bigger_than(otmp);
	
	while((p) and (p->position <= right))
	{
// 		otmp.position = p->position;
		origins.remove(p);
		p = origins.find_bigger_than(otmp);
		cntr++;
	}
}

void replicator::checkall_annihilations() const
{
	unsigned int n_left = rf_left.get_n_elements();
	unsigned int n_right = rf_right.get_n_elements();
	
	for (unsigned int i = 0; i < n_left; i++)
	{
		red_black_tree<fork_t>::const_iterator tmp = rf_left.find_smallest_const(i);
		if (not (annihilations.exists(tmp->anni)))
		{
			throw RuntimeError(string("Annihilation nr ") + tostr(tmp->anni) + string(" does not exist in step " + tostr(step_counter) + string(".")));
		}
	}
	
	for (unsigned int i = 0; i < n_right; i++)
	{
		red_black_tree<fork_t>::const_iterator tmp = rf_right.find_smallest_const(i);
		if (not (annihilations.exists(tmp->anni)))
		{
			throw RuntimeError(string("Annihilation nr ") + tostr(tmp->anni) + string(" does not exist in step " + tostr(step_counter) + string(".")));
		}
	}
}

void replicator::pre_output()
{
	if (runcounter == 0)
	{
		log_file.open(parameters.log_file_name.c_str());
	}
	if (parameters.write_to_disk)
	{
		rate_file.open(parameters.rate_file_name.c_str());
// 		anni_file.open(parameters.anni_file_name.c_str());
// 		lifetime_file.open(parameters.lifetimes_file_name.c_str());
		if(parameters.insane_origin_debug_level)
			ori_d_file.open(parameters.ori_d_file_name.c_str());
		
// 		ofstream ori_file(parameters.ori_file_name.c_str());
// 		for (unsigned int i = 0; i < origins.get_n_elements(); i++)
// 		{
// 			ori_file << std::scientific << (origins.find_biggest(i))->position << " " << (origins.find_biggest(i))->chromatin << endl;
// 		}
// 		ori_file.close();
		
		if (parameters.write_covered_origins)
		{
// 			co_file.open(parameters.covered_origin_file_name.c_str());
		}
		if (parameters.write_induced_patches)
		{
// 			ip_file.open(parameters.ip_file_name.c_str());
		}
		if (parameters.dump_fork_positions)
		{
// 			fork_file.open(parameters.fork_file_name.c_str());
//			char fileName[] = "dada.dat";
			string binaryString = parameters.fork_file_name+"Binary";
			fork_file_binary = fopen(binaryString.c_str(), "wb");
		}
	}
	
	log_file << "Started run Nr. " << runcounter << "." << endl;
}

void replicator::output(double out_time)
{
	fork_consistency_check();
	if (parameters.write_to_disk)
	{
		vector<unsigned int> fcounter;
		for (unsigned int j = 0; j < chromatins.size(); j++)
		{
			fcounter.push_back(0);
		}
		
		unsigned int n_elems = rf_left.get_n_elements();
		for(unsigned int i = 0; i < n_elems; i++)
		{
			red_black_tree<fork_t>::iterator p = rf_left.find_biggest(i);
			
			fcounter[p->chromatin]++;
		}
		
		if (parameters.write_induced_patches)
		{
			// We have to determine the number of induced firing patches, their average size and the average number of origins they cover.
			unsigned int n_ip = 0;
			unsigned int li = 0;
			unsigned int ri = 0;
			unsigned int n_covered = 0;
			unsigned int n_cforks = 0;
			double avg_size = 0.;
			double connection_dist = 1.e6;
			bool curr_left = true;
			bool next_left = false;
			double lo_end = 0.;
			double hi_end = 0.;
			
			red_black_tree<fork_t>::iterator pl = rf_left.find_smallest();
			red_black_tree<fork_t>::iterator pr = rf_right.find_smallest();
			red_black_tree<fork_t>::iterator current = pl;
			red_black_tree<fork_t>::iterator next = pr;
			
			if ((pl != NULL) and (pr != NULL))
			{
				if (*next < *current)
				{
					next = pl;
					current = pr;
					curr_left = false;
					next_left = true;
				}
				
				lo_end = current->position;
				hi_end = lo_end;
			}
			
			
			while ((pl != NULL) and (pr != NULL))
			{
				if (next->position - current->position < connection_dist)
				{
					// Still in the same patch.
					if (curr_left or (not next_left))
					{
						// Count the origins in between them.
						n_covered += count_origins(current->position, next->position);
						n_cforks += count_forks(current->position, next->position);;
					}
				}
				else
				{
					n_ip++;
					avg_size += hi_end-lo_end;
					lo_end = next->position;
					
					if (curr_left or (not next_left))
					{
						// Count the origins in the two covered intervals.
						n_covered += count_origins(current->position, current->position + connection_dist/2.);
						n_covered += count_origins(next->position - connection_dist/2., next->position);
						
						n_cforks += count_forks(current->position, current->position + connection_dist/2.);
						n_cforks += count_forks(next->position - connection_dist/2., next->position);
					}
				}
				hi_end = next->position;
				
				while (pl->position <= next->position)
				{
					li++;
					pl = rf_left.find_smallest(li);
					if (pl == NULL)
						break;
				}
				while (pr->position <= next->position)
				{
					ri++;
					pr = rf_right.find_smallest(ri);
					if (pr == NULL)
						break;
				}
				if ((pl == NULL) or (pr == NULL))
					continue;
				
				current = next;
				curr_left = next_left;
				if (*pl < *pr)
				{
					next = pl;
					next_left = true;
				}
				else
				{
					next = pr;
					next_left = false;
				}
			}
			pl = rf_left.find_smallest();
			pr = rf_right.find_smallest();
			red_black_tree<fork_t>::iterator pl2 = rf_left.find_biggest();
			red_black_tree<fork_t>::iterator pr2 = rf_right.find_biggest();
			
			if (not (pl == NULL or pr == NULL or pl2 == NULL or pr2 == NULL))
			{
				if (min(pl->position, pr->position) - max(pl2->position, pr2->position)+parameters.strand_length > 1.e6)
				{
					n_ip++;
				}
				avg_size += max(pl2->position, pr2->position)-lo_end;
			}
			
			double dn_covered = 0.;
// 			double dn_cforks = 0.;
			if (n_ip == 0)
			{
				avg_size = 0.;
			}
			else
			{
				avg_size /= n_ip;
				dn_covered = double(n_covered)/n_ip;
// 				dn_cforks = double(n_cforks)/n_ip;
			}
			
// 			ip_file << current_time << "\t" << n_ip << "\t" << avg_size << "\t" << dn_covered << "\t" << n_cforks << endl;
		}
		
		n_elems = rf_right.get_n_elements();
		for(unsigned int i = 0; i < n_elems; i++)
		{
			red_black_tree<fork_t>::iterator p = rf_right.find_biggest(i);
			boundary btemp(p->position, 0, 0, 0, 0, 0);
			red_black_tree<boundary>::iterator q = boundaries.find_bigger_than(btemp);
			if (q)
			{
				fcounter[q->left]++;
			}
			else
			{
				red_black_tree<boundary>::iterator q2 = boundaries.find_smallest(0);
				if (q2)
				{
					fcounter[q2->left]++;
				}
			}
			
		}
		
		rate_file << std::scientific << out_time;
		
		for (vector<unsigned int>::iterator pp = fcounter.begin(); pp != fcounter.end(); ++pp)
		{
			rate_file << "\t" << *pp;
		}
		rate_file << "\t" << origins.get_n_elements();
		rate_file << "\t" << n_fired_origins;
		
		for (vector<unsigned int>::iterator pt = trig_firing.begin(); pt != trig_firing.end(); ++pt)
		{
			rate_file << "\t" << *pt;
		}
		
		for (vector<unsigned int>::iterator pu = untrig_firing.begin(); pu != untrig_firing.end(); ++pu)
		{
			rate_file << "\t" << *pu;
		}
		
		for (vector<unsigned int>::iterator pc = n_crossovers.begin(); pc != n_crossovers.end(); ++pc)
		{
			rate_file << "\t" << *pc;
		}
		
		rate_file << "\n";
		
// 		anni_file << std::scientific << current_time;
// 		for (vector<unsigned int>::iterator pa = anni_count.begin(); pa != anni_count.end(); ++pa)
// 		{
// 			anni_file << std::scientific << "\t" << *pa;
// 		}
// 		anni_file << endl;
		
// 		lifetime_file << current_time;
// 		double curtime = 0.;
// 		for(unsigned int j = 0; j < lifetime_count->size(); j++, curtime += parameters.lifetime_bin_width)
// 		{
// 			lifetime_file << "\t" << curtime;
// 			for (vector<unsigned int>::iterator p = (*lifetime_count)[j].begin(); p != (*lifetime_count)[j].end(); ++p)
// 			{
// 				lifetime_file << "\t" << *p;
// 			}
// 		}
// 		
// 		lifetime_file << endl;
		
		if (parameters.write_covered_origins)
		{
// 			co_file << std::scientific << out_time;
// 			pair<vector<long double>, vector<unsigned int> > resp = get_covered_origins();
// 			for (vector<long double>::const_iterator p = resp.first.begin(); p != resp.first.end(); ++p)
// 			{
// 				co_file << std::scientific <<  "\t" << *p;
// 			}
// 			for (vector<unsigned int>::const_iterator q = resp.second.begin(); q != resp.second.end(); ++q)
// 			{
// 				co_file << std::scientific << "\t" << *q;
// 			}
// 			co_file << "\n";
		}
		
		if (parameters.dump_fork_positions)
		{
			// nicor binary
			double f = -1.;
			fwrite(&f, sizeof(double), 1, fork_file_binary);
			// rocin
// 			fork_file << std::scientific << std::setprecision(std::numeric_limits<double>::digits10) << out_time;

			unsigned int n_elems = rf_left.get_n_elements();
			for(unsigned int i = 0; i < n_elems; i++)
			{
				red_black_tree<fork_t>::iterator p = rf_left.find_biggest(i);
// 				fork_file << "\t" << p->position + (out_time - p->time)*chromatins[p->chromatin].fork_speed;
				// nicor binary
 				f = (double) (p->position + (out_time - p->time)*chromatins[p->chromatin].fork_speed);
 				fwrite(&f, sizeof(double), 1, fork_file_binary);
				// rocin
			}
			n_elems = rf_right.get_n_elements();
			for(unsigned int i = 0; i < n_elems; i++)
			{
				red_black_tree<fork_t>::iterator p = rf_right.find_biggest(i);
// 				fork_file << "\t" << p->position - (out_time - p->time)*chromatins[p->chromatin].fork_speed;
				// nicor binary
 				f = (double) (p->position - (out_time - p->time)*chromatins[p->chromatin].fork_speed);
 				fwrite(&f, sizeof(double), 1, fork_file_binary);
				// rocin
			}
// 			fork_file << endl;
		}
		
		if (parameters.write_note_when_covered)
		{
			vector<set<double>::iterator> toremove;
			for (set<double>::iterator p = parameters.note_when_covered.begin(); p != parameters.note_when_covered.end(); ++p)
			{
				fork_t tmpfork(*p, 0, 0, 0, 0, 0);
				boundary tmpbound(*p, false, 0, 0, 0, 0);
				// nicor
				origin tmporigin(*p, 0, 0);
				// rocin
				red_black_tree<fork_t>::iterator ltol = rf_left.find_smaller_than(tmpfork);
				red_black_tree<fork_t>::iterator ltor = rf_right.find_smaller_than(tmpfork);
				red_black_tree<fork_t>::iterator rtol = rf_left.find_bigger_than(tmpfork);
				red_black_tree<fork_t>::iterator rtor = rf_right.find_bigger_than(tmpfork);
				// nicor
				red_black_tree<origin>::iterator rori = origins.find_bigger_than(tmporigin);
				red_black_tree<origin>::iterator lori = origins.find_smaller_than(tmporigin);
//				std::cout << origins.get_n_elements() << endl;
				
				// hier ????? !!!!! abfragen, ob der left/right origin innerhalb des chromosomes ist. falls ja, kann es noch nicht fertig repliziert sein.
//				std::cout << rori->position - *p << endl;
//				std::cout << *p - lori->position << endl;
				// rocin



				bool leftvalid = true;
				bool rightvalid = true;
				// nicor
//				bool chromosomeFinished = false;
				int tmpChromosome=-1;
				// rocin
				
				double ltolpos = -1.;
				double ltorpos = -1.;
				double rtolpos = -1.;
				double rtorpos = -1.;
				
				if (ltol)
				{
					ltolpos = ltol->position + (current_time - ltol->time)*chromatins[ltol->chromatin].fork_speed;
					if (ltolpos > *p)
					{
						leftvalid = false;
						rightvalid = false;
					}
				}
				else
				{
					leftvalid = false;
				}
				
				if (ltor)
					ltorpos = ltor->position - (current_time - ltor->time)*chromatins[ltor->chromatin].fork_speed;
				if (rtol)
					rtolpos = rtol->position + (current_time - rtol->time)*chromatins[rtol->chromatin].fork_speed;
				if (rtor)
				{
					rtorpos = rtor->position - (current_time - rtor->time)*chromatins[rtor->chromatin].fork_speed;
					if (rtorpos < *p)
					{
						leftvalid = false;
						rightvalid = false;
					}
				}
				else
				{
					rightvalid = false;
				}
				
				red_black_tree<boundary>::iterator lbound = boundaries.find_smaller_than(tmpbound);
				if (not lbound)
					lbound = boundaries.find_smallest();
				
				while(not lbound->is_chromosome_boundary)
					lbound = boundaries.find_smaller_than(*lbound);
				
				red_black_tree<boundary>::iterator rbound = boundaries.find_bigger_than(tmpbound);
				if (not rbound)
					rbound = boundaries.find_biggest();
				
				while(not rbound->is_chromosome_boundary)
					rbound = boundaries.find_bigger_than(*rbound);
				
				// nicor, entweder gibt es keine origins mehr, oder es gibt nur noch einen rechten/linken und dieser ist draussen, oder beide sind draussen (hier geht es um origins nicht forks!)
//				if (origins.get_n_elements() == 0 or (!lori and (rbound->position < rori->position)) or (!rori and (lbound->position > lori->position)) or ((rbound->position < rori->position) and (lbound->position > lori->position)))
//					chromosomeFinished = true;

				if (lori and lbound->position < lori->position)
					tmpChromosome = lori->chromosome;
				if (rori and rbound->position > rori->position)
					tmpChromosome = rori->chromosome;
				// rocin
				if (ltol)
				{
					if (lbound->position > ltolpos)
					{
						leftvalid = false;
					}
					else
					{
						if (ltor)
							if (ltolpos < ltorpos)
							{
								leftvalid = false;
							}
					}
				}
				
				if (rtor)
				{
					if (rbound->position < rtorpos)
					{
						rightvalid = false;
					}
					else
					{
						if (rtol)
							if (rtorpos > rtolpos)
							{
								rightvalid = false;
							}
					}
				}
//				std::stringstream ss;
//				for(size_t bb = 0; bb < chromosomeStarted.size(); ++bb)
//				{
//				  if(bb != 0)
//				    ss << ",";
//				  ss << chromosomeStarted[bb];
//				}
//				std::string s = ss.str();
//				std::cout << s << endl;
				// nicor, entweder es gibt noch eine richtige fork links/rechts oder es wurde noch nicht gestartet
				if (leftvalid or rightvalid or (tmpChromosome!=-1 and !chromosomeStarted[tmpChromosome]))	 // or current_time<5.e3
				// rocin
					continue;
				toremove.push_back(p);
			}
//			std::cout << current_time << endl;

			for(vector<set<double>::iterator>::iterator q = toremove.begin(); q != toremove.end(); ++q)
			{
//				std::cout << current_time << endl;
				covered_at_time.push_back(pair<double,double>(**q, current_time));	// current_time
				parameters.note_when_covered.erase(*q);
			}
		}
	}
}

void replicator::post_output()
{
	log_file << "Finished run Nr. " << runcounter << "." << endl;
	log_file.flush();
	
	if (parameters.write_to_disk)
	{
		output(current_time);
		rate_file.close();
// 		anni_file.close();
// 		lifetime_file.close();
		if (parameters.insane_origin_debug_level)
			ori_d_file.close();
		
		if (parameters.write_covered_origins)
		{
// 			co_file.close();
		}
		if (parameters.write_induced_patches)
		{
// 			ip_file.close();
		}
		if (parameters.dump_fork_positions)
		{
// 			fork_file.close();
			fclose(fork_file_binary);
		}
		
		if (parameters.write_fdd)
		{
			fdd_file.open(parameters.fdd_file_name.c_str());
			
			map<unsigned int, unsigned int> firing_dd;
			
			double last = -1.;
			for(set<double>::iterator p = firing_positions.begin(); p != firing_positions.end(); ++p)
			{
				if (last >= 0.)
				{
					double dist = *p-last;
					unsigned int index = dist/parameters.firing_distance_bin_width;
					if (firing_dd.find(index) == firing_dd.end())
					{
						firing_dd[index] = 1;
					}
					else
					{
						firing_dd[index]++;
					}
				}
				last = *p;
			}
			
			for(map<unsigned int, unsigned int>::iterator p = firing_dd.begin(); p != firing_dd.end(); ++p)
			{
				fdd_file << p->first*parameters.firing_distance_bin_width << "\t" << p->second << endl;
			}
			
			fdd_file.close();
		}
		
		ofstream covered_file(parameters.note_when_covered_file_name.c_str());
		for(list<pair<double,double> >::iterator p = covered_at_time.begin(); p != covered_at_time.end(); ++p)
		{
			covered_file << p->first << "\t" << p->second << endl;
		}
		covered_file.close();
	}
}

unsigned int replicator::count_origins(double start, double end)
{
	// Count the number of origins that lie between start and end.
	
	unsigned int n_found = 0;
	origin tmp(start,0.,-1.);
	red_black_tree<origin>::iterator found = origins.find_bigger_than(tmp);
	while(found != NULL)
	{
		if (found->position >= end)
			break;
		
		n_found++;
		found = origins.find_bigger_than(*found);
	}
	return n_found;
}

unsigned int replicator::count_forks(double start, double end)
{
	// Count the number of origins that lie between start and end.
	
	unsigned int n_found = 0;
	fork_t tmp(start, 0., 0, 0, 0, 0);
	red_black_tree<fork_t>::iterator found = rf_left.find_bigger_than(tmp);
	while(found != NULL)
	{
		if (found->position >= end)
			break;
		
		n_found++;
		found = rf_left.find_bigger_than(*found);
	}
	
	found = rf_right.find_bigger_than(tmp);
	while(found != NULL)
	{
		if (found->position >= end)
			break;
		
		n_found++;
		found = rf_right.find_bigger_than(*found);
	}
	
	return n_found;
}

pair<vector<long double>, vector<unsigned int> > replicator::get_covered_origins() const
{
	// Sum up the induced firing probabilities for all origins that are within induced firing range.
	
	vector<long double> covo(chromatins.size());
	vector<unsigned int> orico(2*(chromatins.size()));
	for(size_t i = 0; i < origins.get_n_elements(); i++)
	{
		red_black_tree<origin>::const_iterator p = origins.find_biggest_const(i);
		if (not (check_skip(p)))
		{
			double probval = origin_prob(p).first;
			if (probval > chromatins[p->chromatin].base)
			{
				covo[p->chromatin] += probval;
				orico[p->chromatin]++;
			}
			orico[chromatins.size() + p->chromatin]++;
		}
	}
	return pair<vector<long double>, vector<unsigned int> > (covo, orico);
}

void replicator::fill_origins()
{
	// Fill the origin tree with origins at random positions.
	
	for(unsigned int i = 0; i < parameters.n_origins; i++)
	{
		bool found = false;
		while (not found)
		{
			double tempos = parameters.strand_length*(*real_it++);
			
			size_t tempchrom = 0;
			int temp_chromosome = -1;
			
			boundary btemp(tempos, 0, 0, 0, 0, 0);
			red_black_tree<boundary>::iterator q = boundaries.find_bigger_than(btemp);
			if (q)
			{
				tempchrom = q->left;
				temp_chromosome = q->left_chromosome;
			}
			else
			{
				red_black_tree<boundary>::iterator q2 = boundaries.find_smallest(0);
				if (q2)
				{
					tempchrom = q2->left;
					temp_chromosome = q2->left_chromosome;
				}
			}
			origin tmpori(tempos, tempchrom, temp_chromosome);
			if (origins.find_same(tmpori))
			{
				continue;
			}
			
			origins.insert(tmpori);
			found = true;
		}
	}
}

void replicator::calc_firing_by_limiter_growth()
{
	unsigned int counter = 0;
	while (counter < parameters.capacity)
	{
		double tmp = double(counter) + 0.5;
		double timeval = log(parameters.capacity/(parameters.capacity-tmp))/parameters.rate;
		new_pairs.insert_value(timeval, counter);
		counter++;
	}
}


void replicator::push_left_anni(fork_t &infork, bool preserve_id = false)
{
	double lhrpos(0.);
	
	red_black_tree<fork_t>::iterator lhit_r = rf_right.find_bigger_than(infork);
	
	if (lhit_r)
		lhrpos = lhit_r->position - chromatins[lhit_r->chromatin].fork_speed * (current_time - lhit_r->time);
	
	boundary tbnd(infork.position, false, 0, 0, 0, 0);
	red_black_tree<boundary>::iterator lhit_b = boundaries.find_bigger_than(tbnd);
	
	anni_type atl = FORK_FORK;
	
	double ldistance = 0.;
	double lspeed = 1.e6;
	
	if (lhit_r)
	{
		if (lhit_b)
		{
			infork.chromatin = lhit_b->left;
			if (lhit_b->position < lhrpos)
			{
				ldistance = lhit_b->position - infork.position;
				
				if (not lhit_b->is_chromosome_boundary)
					atl = FORK_BOUNDARY;
				else
					atl = SOLO_REMOVAL;
			}
			else
			{
				ldistance = lhrpos - infork.position;
			}
			lspeed = chromatins[lhit_b->left].fork_speed;
		}
		else
		{
			ldistance = lhrpos - infork.position;
			red_black_tree<boundary>::iterator lhit_t = boundaries.find_smaller_than(tbnd);
			if(lhit_t)
			{
				infork.chromatin = lhit_t->right;
				lspeed = chromatins[lhit_t->right].fork_speed;
			}
			else
			{
				throw RuntimeError("Left fork couldn\'t find any chromatin boundaries.");
			}
		}
	}
	else
	{
		if (lhit_b)
		{
			ldistance = lhit_b->position - infork.position;
			lspeed = chromatins[lhit_b->left].fork_speed;
			infork.chromatin = lhit_b->left;
			
			if (not lhit_b->is_chromosome_boundary)
				atl = FORK_BOUNDARY;
			else
				atl = SOLO_REMOVAL;
		}
		else
		{
			// We have to look at the lower end for a hit.
			fork_t t2f(0., current_time, 0, 0, 0, 0);
			boundary t2b(t2f.position, false, 0, 0, 0, 0);
			lhit_r = rf_right.find_bigger_than(t2f);
			if (lhit_r)
				lhrpos = lhit_r->position - chromatins[lhit_r->chromatin].fork_speed * (current_time - lhit_r->time);
			
			lhit_b = boundaries.find_bigger_than(t2b);
			
			if (not lhit_r)
			{
				if (lhit_b)
				{
					infork.chromatin = lhit_b->left;
					ldistance = lhit_b->position + (parameters.strand_length - infork.position);
					lspeed = chromatins[lhit_b->left].fork_speed;
					
					if (not lhit_b->is_chromosome_boundary)
						atl = FORK_BOUNDARY;
					else
						atl = SOLO_REMOVAL;
				}
				else
				{
					throw RuntimeError("Left fork couldn\'t find any chromatin boundaries.");
				}
			}
			else
			{
				if (lhit_b)
				{
					infork.chromatin = lhit_b->left;
					if (lhit_b->position < lhrpos)
					{
						ldistance = lhit_b->position + (parameters.strand_length - infork.position);
						if (not lhit_b->is_chromosome_boundary)
							atl = FORK_BOUNDARY;
						else
							atl = SOLO_REMOVAL;
					}
					else
					{
						ldistance = lhrpos + (parameters.strand_length - infork.position);
					}
					lspeed = chromatins[lhit_b->left].fork_speed;
				}
				else
				{
					throw RuntimeError("Left fork couldn\'t find any chromatin boundaries.");
				}
			}
		}
	}
	
	// Now build the actual annihilation event to insert.
	annihilation lanni;
	lanni.at = atl;
	if (atl == FORK_FORK)
		lanni.endtime = ldistance/2./lspeed + current_time;
	else
		lanni.endtime = ldistance/lspeed + current_time;
	
	if (not preserve_id)
	{
		infork.id = fork_counter;
		fork_counter++;
	}
	
	lanni.first_fork = rf_left.insert(infork);
	lanni.left = true;
	
	
	if (atl == FORK_FORK)
	{
		lanni.second_fork = lhit_r;
	}
	else
	{
		lanni.other_boundary = lhit_b;
	}
	
	annihilations.insert_value(lanni, anni_counter);
	
	lanni.first_fork->anni = anni_counter;
	if (lanni.at == FORK_FORK)
	{
		if (annihilations.exists(lanni.second_fork->anni))
		{
			annihilations.remove(lanni.second_fork->anni);
		}
		lanni.second_fork->anni = anni_counter;
	}
	anni_counter++;
	
}

void replicator::push_right_anni(fork_t &infork, bool preserve_id = false)
{
	double rhlpos(0.);
	
	red_black_tree<fork_t>::iterator rhit_l = rf_left.find_smaller_than(infork);
	
	if (rhit_l)
		rhlpos = rhit_l->position + chromatins[rhit_l->chromatin].fork_speed * (current_time - rhit_l->time);
	
	boundary tbnd(infork.position, false, 0, 0, 0, 0);
	red_black_tree<boundary>::iterator rhit_b = boundaries.find_smaller_than(tbnd);
	
	anni_type atr = FORK_FORK;
	
	double rdistance = 0.;
	double rspeed = 1.e6;
	
	if (rhit_l)
	{
		if (rhit_b)
		{
			infork.chromatin = rhit_b->right;
			if (rhit_b->position > rhlpos)
			{
				rdistance = infork.position - rhit_b->position;
				
				if (not rhit_b->is_chromosome_boundary)
					atr = FORK_BOUNDARY;
				else
					atr = SOLO_REMOVAL;
			}
			else
			{
				rdistance = infork.position - rhlpos;
			}
			rspeed = chromatins[rhit_b->right].fork_speed;
		}
		else
		{
			rdistance = infork.position - rhlpos;
			red_black_tree<boundary>::iterator rhit_t = boundaries.find_bigger_than(tbnd);
			if(rhit_t)
			{
				infork.chromatin = rhit_t->left;
				rspeed = chromatins[rhit_t->left].fork_speed;
			}
			else
			{
				throw RuntimeError("Right fork couldn\'t find any chromatin boundaries.");
			}
		}
	}
	else
	{
		if (rhit_b)
		{
			infork.chromatin = rhit_b->right;
			rdistance = infork.position - rhit_b->position;
			rspeed = chromatins[rhit_b->right].fork_speed;
			
			if (not rhit_b->is_chromosome_boundary)
				atr = FORK_BOUNDARY;
			else
				atr = SOLO_REMOVAL;
		}
		else
		{
			// We have to look at the lower end for a hit.
			fork_t t2f(parameters.strand_length, current_time, 0, 0, 0, 0);
			boundary t2b(parameters.strand_length, false, 0, 0, 0, 0);
			rhit_l = rf_left.find_smaller_than(t2f);
			if (rhit_l)
				rhlpos = rhit_l->position + chromatins[rhit_l->chromatin].fork_speed * (current_time - rhit_l->time);
			
			rhit_b = boundaries.find_smaller_than(t2b);
			
			if (not rhit_l)
			{
				if (rhit_b)
				{
					infork.chromatin = rhit_b->right;
					rdistance = infork.position + (parameters.strand_length - rhit_b->position);
					rspeed = chromatins[rhit_b->right].fork_speed;
					
					if (not rhit_b->is_chromosome_boundary)
						atr = FORK_BOUNDARY;
					else
						atr = SOLO_REMOVAL;
				}
				else
				{
					throw RuntimeError("Right fork couldn\'t find any chromatin boundaries.");
				}
			}
			else
			{
				if (rhit_b)
				{
					infork.chromatin = rhit_b->right;
					if (rhit_b->position > rhlpos)
					{
						rdistance = infork.position + (parameters.strand_length - rhit_b->position);
						
						if (not rhit_b->is_chromosome_boundary)
							atr = FORK_BOUNDARY;
						else
							atr = SOLO_REMOVAL;
					}
					else
					{
						rdistance = infork.position + (parameters.strand_length - rhlpos);
					}
					rspeed = chromatins[rhit_b->right].fork_speed;
				}
				else
				{
					throw RuntimeError("Left fork couldn\'t find any chromatin boundaries.");
				}
			}
		}
	}
	
	// Now build the actual annihilation event to insert.
	annihilation ranni;
	
	if (not preserve_id)
	{
		infork.id = fork_counter;
		fork_counter++;
	}
	
	ranni.first_fork = rf_right.insert(infork);

	if (atr == FORK_FORK)
		ranni.endtime = rdistance/2./rspeed + current_time;
	else
		ranni.endtime = rdistance/rspeed + current_time;
	
	ranni.at = atr;
	ranni.left = false;
	
	if (atr == FORK_FORK)
	{
		ranni.second_fork = rhit_l;
	}
	else
	{
		ranni.other_boundary = rhit_b;
	}
	annihilations.insert_value(ranni, anni_counter);
		
	ranni.first_fork->anni = anni_counter;
	if (ranni.at == FORK_FORK)
	{
		// Existence check has to be there, since else it is possible to remove an annihilation twice, because it is listed by both of the forks that are changed in the FORK_FORK case.
		if (annihilations.exists(ranni.second_fork->anni))
		{
			annihilations.remove(ranni.second_fork->anni);
		}
		ranni.second_fork->anni = anni_counter;
	}
	anni_counter++;
}

void replicator::fire_origin(red_black_tree<origin>::iterator tofire, unsigned int inid)
{
	if (not tofire)
		return;
	
	firing_positions.insert(tofire->position);
	
	fork_t tmpf(tofire->position, current_time, inid, fork_counter, 0, tofire->chromosome);
	
	push_left_anni(tmpf);
	push_right_anni(tmpf);
	
	// hier ?????

	chromosomeStarted[tofire->chromosome] = true;
	// ?????

	// At the very end, remove the origin that was just fired.
	origins.remove(tofire);
	n_fired_origins++;
}


bool replicator::test_annihilation_consistency()
{
	unsigned int n_annis = annihilations.get_n_elements();
	unsigned int n_forks = rf_left.get_n_elements() + rf_right.get_n_elements();
	
	if (2*n_annis < n_forks)
		return false;
	if (n_annis > n_forks)
		return false;
	return true;
}

void replicator::set_parameters(const parameter_set &params)
{
	parameters = params;
	chromosomeStarted = params.chromosomeStarted;
	max_induced_firing = params.max_induced_firing;
	boundaries.clear();
	if (not (lifetime_count))
		lifetime_count = new vector< vector<unsigned int> >;
	
	bool foundfirst = false;
	bool foundlast = false;
	
	for(vector<boundary>::const_iterator p = params.boundaries.begin(); p != params.boundaries.end(); ++p)
	{
		if (boundaries.find_same(*p))
			continue;
		
		boundaries.insert(*p);
		if ((p->position == parameters.strand_length) or (*p == params.boundaries.back()))
		{
			if (p->is_chromosome_boundary)
				foundlast = true;
		}
		if (p->position == 0)
		{
			if (p->is_chromosome_boundary)
				foundfirst = true;
		}
	}
	
	if (not (foundlast and foundfirst))
		throw RuntimeError("Chromatin range does not exactly start/end with a chromosome boundary.");
	
	chromatins.clear();
	chromatins = params.chromatins;
	
	parameter_changes.clear();
	size_t cntr = 0;
	for(vector<change_parameter>::const_iterator cp = params.change_params.begin(); cp != params.change_params.end(); ++cp, ++cntr)
	{
		parameter_changes.insert_value(*cp,cntr);
	}
	
	for (size_t i = 0; i < chromatins.size(); i++)
	{
		trig_firing.push_back(0);
		untrig_firing.push_back(0);
		anni_count.push_back(0);
	}
	
	n_crossovers.clear();
	for (size_t j = 0; j < (chromatins.size())*(chromatins.size()); j++)
	{
		n_crossovers.push_back(0);
	}
	
	if (parameters.use_origin_list)
	{
		unsigned int tcnt = 0;
		unsigned ori_counter = 0;
		for(vector<double>::iterator p = parameters.origins.begin(); p != parameters.origins.end(); ++p, ++ori_counter)
		{
			size_t tempchrom = 0;
			int temp_chromosome = -1;
				
			boundary btemp(*p, 0, 0, 0, 0, 0);
			red_black_tree<boundary>::iterator q = boundaries.find_bigger_than(btemp);
			if (q)
			{
				tempchrom = q->left;
				temp_chromosome = q->left_chromosome;
			}
			else
			{
				red_black_tree<boundary>::iterator q2 = boundaries.find_smallest(0);
				if (q2)
				{
					tempchrom = q2->left;
					temp_chromosome = q2->left_chromosome;
				}
			}
			origin tmpori(*p, tempchrom, temp_chromosome);
			if (origins.find_same(tmpori))
			{
				std::cout << "DUPLICATE!!" << std::endl;
				continue;
			}
			
			// Discard origins that are malformed (i.e. are on or outside the numerical limits of our genome).
			// The number of such origins should be at most two per calculation, thus negligible. (Unless, of course,
			// the input data is complete garbage.)
			
			if ((*p <= params.boundaries[0].position) or (*p >= params.boundaries[params.boundaries.size()-1].position))
			{
				continue;
			}
			
			red_black_tree<origin>::iterator rboit = origins.insert(tmpori);
			if (parameters.earlyfire.find(ori_counter) != parameters.earlyfire.end())
			{
				fire_origin(rboit, cascade_counter);
				untrig_firing[tmpori.chromatin]++;
				cascade_counter++;
				n_active_pairs++;
				tcnt++;
			}
		}
	}
	else
	{
		fill_origins();
	}
	
	for (unsigned i = 0; i < origins.get_n_elements(); ++i)
	{
		red_black_tree<origin>::iterator p = origins.find_smallest(i);
		boundary tbound(p->position, false, 0, 1, 0, 1);
		red_black_tree<boundary>::iterator q = boundaries.find_smaller_than(tbound);
		p->chromosome = q->right_chromosome;
	}
}

void replicator::store_lifetime(const double lifetime, const size_t chromatin)
{
	unsigned int index = (unsigned int)(lifetime/parameters.lifetime_bin_width);
	while(index + 1 > lifetime_count->size())
	{
		vector<unsigned int> tmp(parameters.chromatins.size());
		lifetime_count->push_back(tmp);
	}
	(*lifetime_count)[index][chromatin]++;
}

bool replicator::check_skip(red_black_tree<origin>::const_iterator p) const
{
	bool skip = false;
	fork_t tmpf(p->position, 0, 0, 0, 0, 0);
		
	red_black_tree<fork_t>::const_iterator lf_hit = rf_left.find_smaller_than_const(tmpf);
	if(lf_hit)
	{
		double lpos = lf_hit->position + chromatins[lf_hit->chromatin].fork_speed*(current_time - lf_hit->time);
		if (lpos >= p->position)
		{
			skip = true;
		}
	}
		
	red_black_tree<fork_t>::const_iterator lf_hit2 = rf_right.find_bigger_than_const(tmpf);
	if (lf_hit2)
	{
		double rpos = lf_hit2->position - chromatins[lf_hit2->chromatin].fork_speed*(current_time - lf_hit2->time);
		if (rpos <= p->position)
		{
			skip = true;
		}
	}
	return skip;
}

pair<double, int> replicator::origin_prob(red_black_tree<origin>::const_iterator p) const
{
	double mindist = 1.e20;
	
	int origin_id = -2; // If there is no fork in this chromosome then return this value. Hier -2 damit ich zwischen "keine fork in der naehe" und "innerhalb inhibition distance" unterscheiden kann
	fork_t tmpf(p->position, 0, 0, 0, 0, 0);
	
	red_black_tree<fork_t>::const_iterator rhit = rf_right.find_bigger_than_const(tmpf);
	red_black_tree<fork_t>::const_iterator lhit = rf_left.find_smaller_than_const(tmpf);
	
	red_black_tree<fork_t>::const_iterator rhit_other = rf_right.find_smaller_than_const(tmpf);
	red_black_tree<fork_t>::const_iterator lhit_other = rf_left.find_bigger_than_const(tmpf);
	
	if (rhit and (rhit->chromosome == p->chromosome))
	{
		double rhlpos = rhit->position - chromatins[rhit->chromatin].fork_speed* (current_time - rhit->time);
		mindist = rhlpos - p->position;
		origin_id = rhit->origin_id;
	}
	
	if (lhit and (lhit->chromosome == p->chromosome))
	{
		double lhrpos = lhit->position + chromatins[lhit->chromatin].fork_speed* (current_time - lhit->time);
		if ((p->position - lhrpos) < mindist)
		{
			mindist = p->position - lhrpos;
			origin_id = lhit->origin_id;
		}
	}
	
	if (rhit_other and (rhit_other->chromosome == p->chromosome))
	{
		double rhlpos = rhit_other->position - chromatins[rhit_other->chromatin].fork_speed* (current_time - rhit_other->time);
		if ((p->position - rhlpos) < mindist)
		{
			mindist = p->position - rhlpos;
			origin_id = rhit_other->origin_id;
		}
	}
	
	if (lhit_other and (lhit_other->chromosome == p->chromosome))
	{
		double lhrpos = lhit_other->position + chromatins[lhit_other->chromatin].fork_speed* (current_time - lhit_other->time);
		if ((lhrpos - p->position) < mindist)
		{
			mindist = lhrpos - p->position;
			origin_id = lhit_other->origin_id;
		}
	}
	
	if (parameters.use_inhibition)
	{
		if (mindist < parameters.inhibition_distance)
			{
				return pair<double, int>(0.,-1);
			}
	}
	
	if (parameters.use_alternate_induced_firing)
	{
		if (mindist < chromatins[p->chromatin].sigma)
			return pair<double, int>(1., origin_id);
		return pair<double, int>(0.,origin_id);
	}

	double gaussian_prob = exp(-mindist*mindist/(2.* chromatins[p->chromatin].sigma* chromatins[p->chromatin].sigma));
	if (parameters.use_induced_limit)
	{
		// In this case there exists a limit for the induced firing probability. If its value lies
		// below that limit, it is set to zero.
		if (chromatins[p->chromatin].induced_limit > gaussian_prob)
		{
			gaussian_prob = 0.;
		}
	}
	
	// Use the standard induced firing function, a cut off gaussian.
	return pair<double, int>(gaussian_prob, origin_id);
// 	return pair<double, int>(exp(-log(mindist/chromatins[p->chromatin].sigma)*log(mindist/chromatins[p->chromatin].sigma)), origin_id);
}

void replicator::fill_by_firing()
{
	while (n_active_pairs < current_limit)
	{
		if (origins.get_n_elements() > 0)
		{
			pair<red_black_tree<origin>::iterator, int> orif = get_firing_origin();
			if (orif.second != -1)	 // orif.second >= 0 hier ????? vielleicht != -1 damit die -2 nicht wegfallen
			{
				fire_origin(orif.first, orif.second);
			}
			else
			{
				// We can't fire any more origins for now.
				break;
			}
		}
		// hier ????? dieses else war im inneren if, aber macht weniger sinn oda?

		n_active_pairs++;
	}
}

pair<red_black_tree<origin>::iterator, int> replicator::get_firing_origin()
{
	// At this point we know that an origin should fire and that there is at least one origin left that can be fired.
	
	bool success = false;
	
	uniform_int_t my_uni(0, origins.get_n_elements() - 1);
	
	// Must use pointer/new construct here, because the assignment made a couple of lines down is 
	// not possible for a uni_int_gen type object. 
	// The reason is that uni_int_gen has a reference member which should under no circumstances
	// be assigned.
	uni_int_gen *ug = new uni_int_gen(rng_engine, my_uni);
	boost::generator_iterator<uni_int_gen> ugi(ug);
	
	unsigned int failure_counter(0);
	while (not success and (origins.get_n_elements() > 0 ))
	{
		if (not firing_prob_summation_mode)
		{
			// Randomly select one of the available origins.
			unsigned int indx = *ugi++;

// 			if (parameters.insane_origin_debug_level)
// 			{
// 				ori_d_file << indx << endl;
// 			}
			
			red_black_tree<origin>::iterator p = origins.find_smallest(indx);
//			cout << std::scientific << indx << endl;
//			cout << std::scientific << p->position << endl;
			
			// Check, if this origin has been swept by a fork. If so, remove it and update the origins container.
			if (check_skip(p))
			{
				origins.remove(p);
				my_uni = uniform_int_t(0, origins.get_n_elements() - 1);
				delete ug;
				ug = new uni_int_gen(rng_engine, my_uni);
				ugi = boost::generator_iterator<uni_int_gen>(ug);
				continue;
			}
			
			pair<double, int> probval = origin_prob(p);
			
			// hier ?????, damit heterochromatin spaeter dran kommt, gauss abschwaechen, dann gibt es aber im euchromatin fast keine induzierten mehr
			// ???????????????? bei 0.3 kommt immer bei 11. simulation segmentation fault, bei 1.0 immer bei der 3. simulation, seltsam!!!, anscheinend beim note_when_covered!!!
			// nicor
			probval.first *= max_induced_firing;
			// rocin
			bool triggered = true;
			if ((probval.second != -1) and (chromatins[p->chromatin].base > probval.first))	// hier ????? warum diese abfrage, dann verhindert man doch sponantes feuern, das hier war noch drin: (probval.second != -1) and
			{
				probval.first = chromatins[p->chromatin].base;
				triggered = false;
			}
			
			// hier ????? wahrscheinlichkeit fuer inaktives x-chromosom
//			if(p->chromosome==75)
//			{
//				probval.first = chromatins[2].base;
//				triggered = false;
//			}

			// nicor: if abfrage nach origin_prob() verschoben
// 			if (parameters.use_induced_limit)
// 			{
// 				// In this case there exists a limit for the induced firing probability. If its value lies
// 				// below that limit, it is set to zero.
// 				if (chromatins[p->chromatin].induced_limit > probval.first)
// 				{
// 					probval.first = 0;
// 				}
// 			}
			// rocin
			
			// Compare probability with a random number to see, if we were successful.
			
			double rn = *real_it++;
			
			// Checking the failure counter is a hack that is needed to avoid an infinite loop in the rare event that all origins that are left have a basal probability of zero and are too far away from all forks.

			if (rn < probval.first)
			{
				delete ug;
				if (parameters.insane_origin_debug_level)
				{
					ori_d_file << std::scientific << current_time << "\t" << p->chromatin << "\t" << triggered << "\t" << p->position << "\t" << failure_counter << endl;
				}
				
				// falls rn unter der induzierten ist, aber auch unter der sponaten w-keit
				// Origins in the induced firing zone can also be fired spontaneously.
				if (rn < chromatins[p->chromatin].base)	// rn < chromatins[p->chromatin].base	// probval.second==-1
				{
					triggered = false;
				}
				
				// Determine, which chromatin type the origin is in.
				
				if (triggered)
				{
					trig_firing[p->chromatin]++;
					return pair<red_black_tree<origin>::iterator, int>(p, probval.second);
				}
				else
				{
					untrig_firing[p->chromatin]++;
					double cctmp = cascade_counter;
					cascade_counter++;
					return pair<red_black_tree<origin>::iterator, int>(p, cctmp);
				}
			}
			else
			{
				failure_counter++;
				if (failure_counter > 1e6)
				{
					firing_prob_summation_mode = true;
					log_file << "Going to fps mode at " << current_time << "." << endl;
					cout << std::scientific << origins.get_n_elements() << " ultra slow mode" << endl;
				}
			}
		}
		else
		{
			// We are in the ultraslow probability summation mode, made necessary by the fact, that the main mode breaks down if all remaining origins have a spontaneous firing probability of zero.

//			if(gaussProb<.2)
//			{
//				gaussProb+= 0.02;
//			}
//			gaussProb = 0.2;
			
			vector<double> temp;
			temp.reserve(origins.get_n_elements());
			
			// While building the probability tree also clean out all the origins because in probability_summation mode there may not be any covered origins.
			unsigned int i = 0;
			while(i < origins.get_n_elements())
			{
				red_black_tree<origin>::iterator q = origins.find_smallest(i);
				
						// Check, if this origin has been swept by a fork. If so, remove it and update the origins container.
				if (check_skip(q))
				{
					origins.remove(q);
					my_uni = uniform_int_t(0, origins.get_n_elements() - 1);
					delete ug;
					ug = new uni_int_gen(rng_engine, my_uni);
					ugi = boost::generator_iterator<uni_int_gen>(ug);
					continue;
				}
				else
				{
					// Add this entries probability value to the list.
					pair<double, int> probval = origin_prob(q);
			
					if ((probval.second != -1) and (chromatins[q->chromatin].base > probval.first))	// hier ?????  (probval.second != -1) and
					{
						probval.first = chromatins[q->chromatin].base;
					}

					// nicor: if abfrage nach origin_prob() verschoben
// 					if (parameters.use_induced_limit)
// 					{
// 						if(chromatins[q->chromatin].induced_limit > probval.first)
// 							probval.first = 0.;
// 					}
					// rocin
					
					temp.push_back(probval.first);
					i++;
				}
			}
			
			if (origins.get_n_elements() == 0)
				continue;
			
			binary_sum_search<double> sums(temp);
			
			double fullprop = sums.get_propensity();
			
// 			if (fullprop < 1.e-7)
			if (fullprop < 0.1)
			{
				return pair<red_black_tree<origin>::iterator, int>(red_black_tree<origin>::iterator(NULL), -1);
			}
			if (fullprop > 1./100.*origins.get_n_elements())
			{
				// If the induced firing probability recovers, kick the simulation
				// back out of ultraslow mode.
				firing_prob_summation_mode = false;
				log_file << "Going back." << endl;
			}
			// Get_index counts from the bottom up, so that we have to use find_smallest.
			unsigned int index = sums.get_index(fullprop*(*real_it++));
			red_black_tree<origin>::iterator toret = origins.find_smallest(index);
			
			double hitval = sums.get_value_from_index(index);
			double rn2 = hitval*(*real_it++);
			
			delete ug;
			if(rn2 <= chromatins[toret->chromatin].base)
			{
				untrig_firing[toret->chromatin]++;
				double cctmp = cascade_counter;
				cascade_counter++;
				return pair<red_black_tree<origin>::iterator, int>(toret, cctmp);
			}
			else
			{
				trig_firing[toret->chromatin]++;
				pair<double, int> probval = origin_prob(toret);
				return pair<red_black_tree<origin>::iterator, int>(toret, probval.second);
			}
		}
		// If no origin was found, repeat.
	}
	delete ug;
	return pair<red_black_tree<origin>::iterator, int>(red_black_tree<origin>::iterator(NULL), -1);
}
