
#include "dna_metropolis.hh"
#include <cmath>

using std::sin;
using std::cos;
using std::exp;
using std::sqrt;
using std::ostream;
using std::endl;
using std::map;
using std::set;
using std::pair;
using std::vector;
using std::list;
using std::cerr;

dna_metropolis::dna_metropolis(unsigned int in_n_beads, double i_kappa_het, double i_kappa_eu, double i_kappa_fac, double i_kappa_2, double i_kappa_3, double i_kappa_fac_rep, double temperature)
  : ellipsis_x(7500), ellipsis_y(5000), ellipsis_z(3500), n_beads(in_n_beads), kappa_het(i_kappa_het), kappa_eu(i_kappa_eu), kappa_fac(i_kappa_fac), kappa_2(i_kappa_2), kappa_3(i_kappa_3), kappa_fac_rep(i_kappa_fac_rep), beta(1.e-18/(boltzmann_constant*temperature)), U(0.), n_cur_step(0), n_accepted(0), use_loops(true), use_gravity(true), start_big_random(false), nonucleus(false), with_nucleoli(true), surface_rad(1000.), hetFactor(0.5), nucleolusRad(1500.), nucleolusRad2(1200.), nucleolusPos({750., 3000., 0.}), nucleolusPos2({-1000., -2500., 0.}), rng_engine(static_cast<std::tr1::mt19937::result_type>(23)), rng_real(0., 1.), vg_real(rng_engine, rng_real), real_it(&vg_real), rng_int(0, n_beads-1), vg_int(rng_engine, rng_int), int_it(&vg_int)
{
}

void dna_metropolis::init(const vector<unsigned> &n_chromosomes, const vector<unsigned> &chromo_beads, const vector<double> &real_ends, const vector<list<pair<size_t, double> > > &band_data, unsigned eu_cons, unsigned het_cons, unsigned fac_cons, unsigned inter_cons, bool truly_random_connections)
{
	srand (time(NULL));
	coordinates.reserve(3*n_beads);
	bead_info.reserve(n_beads);
	
	// The very first thing to do is to generate temporary bead_info values that contain the information on chromosome membership and starters/enders.
	// The connection part of all values is filled with 131071, which is the no_connection indicator.
	unsigned cntr = 0;
	unsigned fc_cntr = 0;
// 	chromosome_starters.insert(0);
	unsigned temppos = 0;
	unsigned oltempos = 0;
	double curx = 0.;
	double oldx = 0.;
	unsigned int no_connection = connection_indicator << 15;
	vector<list<pair<size_t, double> > >::const_iterator bdp = band_data.begin();
	for(vector<unsigned>::const_iterator p = n_chromosomes.begin(); p != n_chromosomes.end(); ++p, ++bdp)
	{
		for (unsigned i = 0; i < *p; ++i)
		{
			list<pair<size_t, double> >::const_iterator lp = bdp->begin();
			oltempos = temppos;
			oldx = curx;
			temppos += chromo_beads[cntr] + 1;
			curx += real_ends[cntr];

			// nicor
// 			bool useStraingingData = true;
// 			unsigned int chromType = 0;
// 			if(usePercent){
// 				double what_type = (double)rand()/(double)RAND_MAX;
// 				if (what_type <= hetPercent)	// heterochromatin wahrscheinlichkeit
// 					chromType = 1;
// 				else if (what_type <= hetPercent + facPercent)	// facultative wahrscheinlichkeit
// 					chromType = 2;
// 			}
			// rocin
			
			pair<unsigned, double> tpx(oltempos, oldx);
			position_offsets.insert(pair<double, pair<unsigned, double> >(curx,tpx));
			for (unsigned j = 0; j < chromo_beads[cntr]+1; ++j)
			{
				if (j*dna_unit > lp->second)
				{
					++lp;
					if (lp == bdp->end())
					{
						--lp;
					}
				}
				unsigned temp_info = fc_cntr;
				if (j == 0)
				{
					temp_info += 128;
				}
				if (j == chromo_beads[cntr])
				{
					temp_info += 256;
				}
				// hier ?????
				// bead_info speichert einfach fuer jedes monomer die information, ob es start/ende ist und welcher chromatin-typ es ist, also keine position oda so.
				// speater wird ueber if((binf>>10)%2) auf heterochromatin gepruft. der >> operator verschiebt die bits nach rechts, die verschobenen fallen raus, also aus 0101 wird 0010, bei 1024 bleibt mit >>10 genau 1 uebrig.
				// Store the chromatin type of this bead in the 10th bit
// 				if(usePercent){
// 				// nicor: falls ich die anteile an eu/het/fac zufaellig bestimmen will
// 					if(chromType==1)
// 					{
// 						temp_info += 1024;
// 					}
// 					// nicor
// 					else if(chromType==2)
// 					{
// 						temp_info += 512;
// 					}
// 				// rocin
// 				} else {
				if(lp->first==1)
				{
					temp_info += 1024;
				}
				// nicor
				else if(lp->first==2)
				{
					temp_info += 512;
				}
// 				}
				
				
				//if (not(j%connection_delta))
				//{
				//	temp_info += 512;
				//	lastloop = j;
				//	if (j == 0)
				//		temp_info += 1024; // Only part of right side loop.
				//	if (chromo_beads[cntr] - j < connection_delta)
				//		temp_info += 2048; // Only part of the left side loop.
				//}
				temp_info += no_connection;
				bead_info.push_back(temp_info);
			}
			fc_cntr++;
		}
		cntr++;
	}
	
	// Now add the n_connections inter-bead connections between beads of the same chromosome. Every bead can have only one connection to another bead and there cannot be connections between beads of different chromosomes. 
	
	unsigned het_cntr = 0;
	unsigned eu_cntr = 0;
	unsigned fac_cntr = 0;
	unsigned inter_cntr = 0;
	
	for (unsigned i = 0; i < het_cons+eu_cons+fac_cons+inter_cons; ++i)
	{
		bool found = false;
		while (not found)
		{
			unsigned k1 = (*int_it++);
			unsigned k2 = (*int_it++);
			
			if (k1 == k2)
				continue;
			if ((bead_info[k1]%128) != (bead_info[k2]%128))
				continue;
			
			// Now we know that they are different beads on the same chromosome.
			// We still have to check, if one of them is already involved in another loop.
			if ((bead_info[k1]>>15 != connection_indicator) or (bead_info[k2]>>15 != connection_indicator))
				continue;
			
			// If we are interested in totally random connections, just accept it.
			if(truly_random_connections)
			{
				// If we are interested in totally random connections, just accept the selected pair.
			}
			else
			{
				// hier ????? facultatives ergaenzen. fac_cons hinzufuegen und: or (bead_info[k1]>>9)%2 != (bead_info[k2]>>9)%2
				// aber braeuchte ich nur ohne truly_random_connections
				// Finally see if there is still capacity for a connection between the two kinds of chromatin.
				if((bead_info[k1]>>10)%2 != (bead_info[k2]>>10)%2 or (bead_info[k1]>>9)%2 != (bead_info[k2]>>9)%2)
				{
					if (inter_cntr < inter_cons)
						++inter_cntr;
					else
						continue;
				}
				else
				{
					if((bead_info[k1]>>10)%2)
					{
						if(het_cntr < het_cons)
							++het_cntr;
						else
							continue;
					}
					else if((bead_info[k1]>>9)%2)
					{
						if(fac_cntr < fac_cons)
							++fac_cntr;
						else
							continue;
					}
					else
					{
						if(eu_cntr < eu_cons)
							++eu_cntr;
						else
							continue;
					}
				}
			}
			
			bead_info[k1] = bead_info[k1]%32768 + (k2 << 15);
			bead_info[k2] = bead_info[k2]%32768 + (k1 << 15);
			found = true;
		}
	}
	
	build_chromatin();
	prepare_chromosomes();
	U = calc_energy();
}

void dna_metropolis::load_coordinates(const vector<double> &in_coords)
{
	coordinates = in_coords;
}

void dna_metropolis::build_chromatin()
{
	// Place beads randomly at fixed distances for the initial chromatin configuration.
	// hier ????? startwert darf nicht im nucleolus liegen! siehe nucleolusPos
	double curx = -500.;
	double cury = 0.;
	double curz = 0.;
	
	for (unsigned i = 0; i < n_beads; ++i)
	{
		coordinates.push_back(curx);
		coordinates.push_back(cury);
		coordinates.push_back(curz);
		
		bool found = false;
		double tcurs[3] = {0., 0., 0.};
		if (start_big_random)
		{
		  	curx = ((*real_it++)*2 -1.)*ellipsis_x;
		  	cury = ((*real_it++)*2 -1.)*ellipsis_y;
		  	curz = ((*real_it++)*2 -1.)*ellipsis_z;
		}
		else
		{
			while (not found)
			{
				double theta = (*real_it++)*pi;
				double psi = (*real_it++)*pi*2;
				double ct = cos(theta);
				double st = sin(theta);
				double cp = cos(psi);
				double sp = sin(psi);
				
				tcurs[0] = curx + st*cp*init_dist;
				tcurs[1] = cury + st*sp*init_dist;
				tcurs[2] = curz + ct*init_dist;
				
				if (isinside(&tcurs[0]))
					found = true;
			}
			curx = tcurs[0];
			cury = tcurs[1];
			curz = tcurs[2];
		}
	}
}

void dna_metropolis::prepare_chromosomes()
{
	if (not use_gravity)
		return;
	
	unsigned cntr = 0;
	
	chromosome tchrom;
	double weight = 0.;
	double xtmp = 0.;
	double ytmp = 0.;
	double ztmp = 0.;
	
	while (cntr < bead_info.size())
	{
		if((bead_info[cntr]>>7)%2)
		{
			// starter
			tchrom.start_pos = cntr;
			tchrom.index = bead_info[cntr]%128;
		}
		xtmp += coordinates[3*cntr];
		ytmp += coordinates[3*cntr+1];
		ztmp += coordinates[3*cntr+2];
		weight += 1.;
		if((bead_info[cntr]>>8)%2)
		{
			// ender
			tchrom.end_pos = cntr;
			chromo_positions.push_back(xtmp/weight);
			chromo_positions.push_back(ytmp/weight);
			chromo_positions.push_back(ztmp/weight);
			chromo_weights.push_back(weight);
			chromosomes.insert(tchrom);
			weight = 0.;
			xtmp = 0.;
			ytmp = 0.;
			ztmp = 0.;
		}
		++cntr;
	}
}

double dna_metropolis::calc_centers_energy() const
{
	double tmpU = 0.;
	
	for(unsigned i = 0; i < ((chromo_weights.size() > 0)?chromo_weights.size()-1:0); ++i)
	{
		for(unsigned j = i+1; j < chromo_weights.size(); ++j)
		{
			double dx = chromo_positions[3*i] - chromo_positions[3*j];
			double dy = chromo_positions[3*i+1] - chromo_positions[3*j+1];
			double dz = chromo_positions[3*i+2] - chromo_positions[3*j+2];
			tmpU += kappa_3*chromo_weights[i]*chromo_weights[j]/(sqrt(dx*dx+dy*dy+dz*dz));
		}
	}
	
	return tmpU;
}

bool dna_metropolis::isinside(const double *p) const
{
	if (nonucleus)
	{
		std::cout << "Ho! " << n_cur_step << std::endl;
		return true;
	}
	// teste, ob punkt in zellkern und ob er ausserhalb des nucleolus ist
	// muss startwert in dna_metropolis::build_chromatin() so legen, dass er nicht im nucleolus liegt
	bool inNucleus = false;
	bool outNucleolus = false;
	
	double xn = *p-nucleolusPos[0];
	double xn2 = *p-nucleolusPos2[0];
	double xv = *p/ellipsis_x;
	++p;
	double yn = *p-nucleolusPos[1];
	double yn2 = *p-nucleolusPos2[1];
	double yv = *p/ellipsis_y;
	++p;
	double zn = *p-nucleolusPos[2];
	double zn2 = *p-nucleolusPos2[2];
	double zv = *p/ellipsis_z;
	
	if (xv*xv+yv*yv+zv*zv <= 1.)
		inNucleus = true;
	
	if (!with_nucleoli)
		return inNucleus;
	if (xn*xn+yn*yn+zn*zn >= nucleolusRad*nucleolusRad && xn2*xn2+yn2*yn2+zn2*zn2 >= nucleolusRad2*nucleolusRad2)
		outNucleolus = true;

	return outNucleolus && inNucleus;

	
}

double dna_metropolis::calc_energy()
{
	// Calculates total energy of system from scratch
	double tmpU = 0.;
	for(unsigned i = 0; i < n_beads-1; ++i)
	{
		if (not ((bead_info[i] >> 8)%2))
		{
			tmpU += en_calc(&coordinates[3*i], &coordinates[3*(i+1)], bead_info[i]);
			tmpU += fac_repuls(&coordinates[3*i], bead_info[i]);	// hier ????? letzte bead nicht dabei, aber die ist eh nicht fakultativ
			if (use_loops)
			{
				unsigned otherind = bead_info[i] >> 15;
				if (otherind != connection_indicator)
				{
					if(otherind > i)
						tmpU += inter_loop_en_calc(&coordinates[3*i], &coordinates[3*(otherind)]);
				}
// 				if (not((bead_info[i]>>11)%2))
// 					tmpU += inter_loop_en_calc(&coordinates[3*i], &coordinates[3*(i+connection_delta)]);
			}
		}
	}
	if (use_gravity)
		tmpU += calc_centers_energy();
	return tmpU;
}

double dna_metropolis::en_calc(const double *left, const double *right, const unsigned &binf) const
{
	// Helper function to calculate energy contribution by bond between left and right bead.
	double dx = *right - *left;
	++right;
	++left;
	double dy = *right - *left;
	++right;
	++left;
	double dz = *right - *left;
	
	if((binf>>10)%2)
		return kappa_het/2.*(dx*dx + dy*dy + dz*dz);
	if((binf>>9)%2)
		return kappa_fac/2.*(dx*dx + dy*dy + dz*dz);
	return kappa_eu/2.*(dx*dx + dy*dy + dz*dz);
}

double dna_metropolis::fac_repuls(const double *left, const unsigned &binf) const
{
	if (!with_nucleoli)
		return 0.;
	
	double xx = *left;
	++left;
	double yy = *left;
	++left;
	double zz = *left;

	
	// fuer alle punkte gilt, dass sie nicht im nucleolus liegen duerfen
	// fuer facultatives gilt zusaetzlich die folgende anziehung
	if((binf>>9)%2)
	{
		// im ersten geht der normierte radius ein, daher ein faktor (ellipsis_x*ellipsis_y*ellipsis_z)**-3
		double membranGrav = -kappa_fac_rep/(surface_rad+5082.*(1.-sqrt(xx/ellipsis_x*xx/ellipsis_x+yy/ellipsis_y*yy/ellipsis_y+zz/ellipsis_z*zz/ellipsis_z)));
		double nucleolusGrav = -kappa_fac_rep/sqrt((xx-nucleolusPos[0])*(xx-nucleolusPos[0])+(yy-nucleolusPos[1])*(yy-nucleolusPos[1])+(zz-nucleolusPos[2])*(zz-nucleolusPos[2]));
		nucleolusGrav += -kappa_fac_rep/sqrt((xx-nucleolusPos2[0])*(xx-nucleolusPos2[0])+(yy-nucleolusPos2[1])*(yy-nucleolusPos2[1])+(zz-nucleolusPos2[2])*(zz-nucleolusPos2[2]));
		return membranGrav + nucleolusGrav;
	}
	// falls ich fuer heterocrhomatin eine abstossung will
	if((binf>>10)%2)
	{
		double membranGrav = hetFactor*kappa_fac_rep/(surface_rad+5082.*(1.-sqrt(xx/ellipsis_x*xx/ellipsis_x+yy/ellipsis_y*yy/ellipsis_y+zz/ellipsis_z*zz/ellipsis_z)));
		double nucleolusGrav = hetFactor*kappa_fac_rep/sqrt((xx-nucleolusPos[0])*(xx-nucleolusPos[0])+(yy-nucleolusPos[1])*(yy-nucleolusPos[1])+(zz-nucleolusPos[2])*(zz-nucleolusPos[2]));
		nucleolusGrav += hetFactor*kappa_fac_rep/sqrt((xx-nucleolusPos2[0])*(xx-nucleolusPos2[0])+(yy-nucleolusPos2[1])*(yy-nucleolusPos2[1])+(zz-nucleolusPos2[2])*(zz-nucleolusPos2[2]));
		return membranGrav + nucleolusGrav;
	}
	return 0.;
}



// double dna_metropolis::en_calc(const double *left, const double *middle, const double *right) const
// {
// 	// Calculate the sum of all energy quantities that this bead contributes to.
// 	double tmpU = 0.;
// 	tmpU += en_calc(left, middle);
// 	tmpU += en_calc(middle, right);
// 	return tmpU;
// }

double dna_metropolis::inter_loop_en_calc(const double *left, const double *right) const
{
	double dx = *right - *left;
	++right;
	++left;
	double dy = *right - *left;
	++right;
	++left;
	double dz = *right - *left;
	
	return kappa_2/2.*(dx*dx+dy*dy+dz*dz);
}

// double dna_metropolis::inter_loop_en_calc(const double *left, const double *middle, const double *right) const
// {
// 	double tmpU = 0.;
// 	tmpU += inter_loop_en_calc(left, middle);
// 	tmpU += inter_loop_en_calc(middle, right);
// 	return tmpU;
// }

double dna_metropolis::center_delta_en_calc(const double *coords, unsigned index, unsigned chromo_index) const
{
	double newx = chromo_positions[3*chromo_index] + (coords[0]-coordinates[3*index])/chromo_weights[chromo_index];
	double newy = chromo_positions[3*chromo_index+1] + (coords[1]-coordinates[3*index+1])/chromo_weights[chromo_index];
	double newz = chromo_positions[3*chromo_index+2] + (coords[2]-coordinates[3*index+2])/chromo_weights[chromo_index];
	
	double delta_U = 0.;
	for(unsigned i = 0; i < chromo_weights.size(); ++i)
	{
		if (i == chromo_index)
			continue;
		
		double dx = chromo_positions[3*i] - newx;
		double dy = chromo_positions[3*i+1] - newy;
		double dz = chromo_positions[3*i+2] - newz;
		
		double dx2 = chromo_positions[3*i] - chromo_positions[3*chromo_index];
		double dy2 = chromo_positions[3*i+1] - chromo_positions[3*chromo_index+1];
		double dz2 = chromo_positions[3*i+2] - chromo_positions[3*chromo_index+2];
		delta_U += kappa_3*chromo_weights[i]*chromo_weights[i]*(1./(sqrt(dx*dx+dy*dy+dz*dz)) - 1./(sqrt(dx2*dx2+dy2*dy2+dz2*dz2)));
		
	}
	return delta_U;
}

double dna_metropolis::get_avg_dist() const
{
	double ret = 0.;
	unsigned count = 0;
	for (unsigned i = 1; i < n_beads; ++i)
	{
		if(bead_info[i-1]%128 != bead_info[i]%128)
			continue;
		double dx = coordinates[3*i] - coordinates[3*(i-1)];
		double dy = coordinates[3*i+1] - coordinates[3*(i-1)+1];
		double dz = coordinates[3*i+2] - coordinates[3*(i-1)+2];
		
		ret += sqrt(dx*dx+dy*dy+dz*dz);
		++count;
	}
	return ret/count;
}

vertex dna_metropolis::dna_to_position(double dna_pos) const
{
	vertex ret;
	
	map<double, pair<unsigned, double> >::const_iterator p = position_offsets.upper_bound(dna_pos);
	
	if (p == position_offsets.end())
	{
		cerr << "Out of boundary position requested in dna_to_position." << endl;
		std::cout << dna_pos << endl;
	}
	
	unsigned delta_pos = unsigned((dna_pos - p->second.second)/dna_unit);
	unsigned pos = p->second.first + delta_pos;
	double valvac = (dna_pos - dna_unit*unsigned(dna_pos/dna_unit))/dna_unit;
	ret.x = coordinates[3*(pos+1)]*valvac + coordinates[3*pos]*(1-valvac);
	ret.y = coordinates[3*(pos+1)+1]*valvac + coordinates[3*pos+1]*(1-valvac);
	ret.z = coordinates[3*(pos+1)+2]*valvac + coordinates[3*pos+2]*(1-valvac);
	return ret;
}

size_t dna_metropolis::dna_to_chromatin_type(double dna_pos) const
{
	map<double, pair<unsigned, double> >::const_iterator p = position_offsets.upper_bound(dna_pos);
	
	if (p == position_offsets.end())
		cerr << "Out of boundary position requested in dna_to_position." << endl;
	
	unsigned delta_pos = unsigned((dna_pos - p->second.second)/dna_unit);
	unsigned pos = p->second.first + delta_pos;
	
	size_t chrom_type = 0;
	if ((bead_info[pos]>>10)%2)
		chrom_type = 1;
	if ((bead_info[pos]>>9)%2)
		chrom_type = 2;
	return chrom_type;
}

double dna_metropolis::run(unsigned n_steps)
{
	// Run n metropolis steps on the system. Return step acceptance rate.
	unsigned n_accepted = 0;
	for(unsigned j = 0; j < n_steps; ++j)
	{
	 	// Regularily re-calculate the total energy to reduce the total error introduced by floating point additions.
		if (not(n_cur_step % 1000000))
		{
			chromo_positions.clear();
			chromo_weights.clear();
			chromosomes.clear();
			prepare_chromosomes();
			U = calc_energy();
		}
		unsigned index = *int_it++;
		
		double newpos[3];
		
		bool foundin = false;
		while (not foundin)
		{
			double theta = (*real_it++)*pi;
			double psi = (*real_it++)*pi*2;
			double ct = cos(theta);
			double st = sin(theta);
			double cp = cos(psi);
			double sp = sin(psi);
			
			newpos[0] = coordinates[3*index] + st*cp*delta_jump;
			newpos[1] = coordinates[3*index+1] + st*sp*delta_jump;
			newpos[2] = coordinates[3*index+2] + ct*delta_jump;
			
			if (isinside(newpos))
				foundin = true;
		}
		
		double U_orig = 0.;
		double U_new = 0.;
		if ((bead_info[index]>>7)%2)
		{
			if ((bead_info[index]>>8)%2)
			{
				// Do nothing since this move makes no difference in energy contribution.
			}
			else
			{
				U_orig = en_calc(&coordinates[3*index], &coordinates[3*(index+1)], bead_info[index]);
				U_orig += fac_repuls(&coordinates[3*index], bead_info[index]);
				U_new = en_calc(&newpos[0], &coordinates[3*(index+1)], bead_info[index]);
				U_new += fac_repuls(&newpos[0], bead_info[index]);
			}
		}
		else if ((bead_info[index]>>8)%2)
		{
			U_orig = en_calc(&coordinates[3*(index-1)], &coordinates[3*index], bead_info[index-1]);
			U_orig += fac_repuls(&coordinates[3*index], bead_info[index]);
			U_new = en_calc(&coordinates[3*(index-1)], &newpos[0], bead_info[index-1]);
			U_new += fac_repuls(&newpos[0], bead_info[index]);
		}
		else
		{
			U_orig = en_calc(&coordinates[3*(index-1)], &coordinates[3*index], bead_info[index-1]);
			U_orig += en_calc(&coordinates[3*index], &coordinates[3*(index+1)], bead_info[index]);
			U_orig += fac_repuls(&coordinates[3*index], bead_info[index]);	// hier ????? abstossung des fakultativen chromatins vom mittelpunkt
			U_new = en_calc(&coordinates[3*(index-1)], &newpos[0], bead_info[index-1]);
			U_new += en_calc(&newpos[0], &coordinates[3*(index+1)], bead_info[index]);
			U_new += fac_repuls(&newpos[0], bead_info[index]);
		}
		
		if (use_loops)
		{
			unsigned tind = bead_info[index]>>15;
			if (tind != connection_indicator)
			{
				U_orig += inter_loop_en_calc(&coordinates[3*tind], &coordinates[3*index]);
				U_new += inter_loop_en_calc(&coordinates[3*tind], &newpos[0]);
			}
		}
		
		double deltaU = U_new - U_orig;
		
		unsigned chromo_index = bead_info[index]%128;
		if (use_gravity)
		{
			deltaU += center_delta_en_calc(&newpos[0], index, chromo_index);
		}
		
		bool accepted = false;
		if (deltaU < 0)
			accepted = true;
		else
		{
			double accprob = exp(-deltaU*beta);
			if ((*real_it++) < accprob)
				accepted = true;
		}
		if (accepted)
		{
			n_accepted++;
			U += deltaU;
			
			if (use_gravity)
			{
				chromo_positions[3*chromo_index] = chromo_positions[3*chromo_index] + (newpos[0]-coordinates[3*index])/chromo_weights[chromo_index];
				chromo_positions[3*chromo_index+1] = chromo_positions[3*chromo_index+1] + (newpos[1]-coordinates[3*index+1])/chromo_weights[chromo_index];
				chromo_positions[3*chromo_index+2] = chromo_positions[3*chromo_index+2] + (newpos[2]-coordinates[3*index+2])/chromo_weights[chromo_index];
			}
			
			coordinates[3*index] = newpos[0];
			coordinates[3*index+1] = newpos[1];
			coordinates[3*index+2] = newpos[2];
		}
		
// 		std::cout << std::scientific << deltaU << "\t" << accepted << endl;
		n_cur_step++;
	}
	return double(n_accepted)/double(n_steps);
}

double dna_metropolis::inter_point_distance(double posa, double posb)
{
	vertex a = dna_to_position(posa);
	vertex b = dna_to_position(posb);
	
	double dist = sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
	return dist;
}

vector<pair<double, double> > dna_metropolis::get_distance_pairs(unsigned n_pairs, double max_dist, double max_position)
{
	vector<pair<double, double> > ret;
	ret.reserve(n_pairs);
	for (unsigned k = 0; k < n_pairs; ++k)
	{
		try
		{
			bool found = false;
			while (!found)
			{
				double delta = (*real_it++)*max_dist;
				double posa = (*real_it++)*max_position;
				double posb = posa + delta;
				
				if (position_offsets.upper_bound(posa)->second.first == position_offsets.upper_bound(posb)->second.first)
				{
					found = true;
					double dist = inter_point_distance(posa,posb);
					ret.push_back(pair<double,double>(delta,dist));
				}
			}
		}
		catch(...)
		{
			cerr << "Exception occurred" << endl;
		}
	}
	return ret;
}

ostream& operator<<(ostream &the_stream, const dna_metropolis &dnam)
{
	unsigned index = 0;
	
	while (index < dnam.n_beads)
	{
		if ((dnam.bead_info[index]>>8)%2)
			the_stream << dnam.coordinates[3*index] << "\t" << dnam.coordinates[3*index+1] << "\t" << dnam.coordinates[3*index+2] << endl;
		else
			the_stream << dnam.coordinates[3*index] << "\t" << dnam.coordinates[3*index+1] << "\t" << dnam.coordinates[3*index+2] << "\t";
		index++;
	}
	return the_stream;
}
