
#include <vector>
#include <iostream>
#include <fstream>
#include "dna_metropolis.hh"
#include <utility>
#include <ctime>
#include <string>
#include <boost/lexical_cast.hpp>

using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::pair;
using std::clock;
using std::string;

int main(void)
{
	vector<unsigned> n_chromosomes = {5,3,3,2,6,3,5,3,5,3,3,3,3,3,3,3,4,2,3,3,3,3,2};
	vector<unsigned> n_beads = {2259, 2426, 1987, 1912, 1795, 1707, 1581, 1462, 1169, 1348, 1328, 1320, 1121, 1047, 983, 875, 786, 754, 614, 620, 446, 488, 1540, 271};
	vector<double> real_ends = {225949719.0, 242651149.0, 198701828.0, 191273063.0, 179557866.0, 170799992.0, 158121424.0, 146274826.0, 116973252.0, 134874737.0, 132852384.0, 132049534.0, 112142981.0, 104768586.0, 98338916.0, 87527254.0, 78674742.0, 75417153.0, 61411651.0, 62035964.0, 44644324.0, 48891433.0, 154013754.0, 27100000.0};
	
	unsigned beadnr = 0;
	unsigned cntr = 0;
	for(vector<unsigned>::iterator p = n_chromosomes.begin(); p != n_chromosomes.end(); ++p)
	{
		beadnr += (*p)*n_beads[cntr] + 1; // Plus one because of the skipped bead in between two chromosomes.
		cntr++;
	}
// 	beadnr--; // Last one needs no skipped bead. // Cancelled out by additional first bead.
	
	dna_metropolis dnam(beadnr, 4.e-7, 4.e-5, 4., 290);
	
	dnam.use_loops = false;
	dnam.use_gravity = false;
	dnam.init(n_chromosomes, n_beads, real_ends);
	
	ofstream tmpfile("tmpfile.txt");
// 		cout << std::scientific << dnam.get_energy() << "\t" << dnam.get_acceptance_ratio() << "\t" << dnam.get_avg_dist() << endl;
	tmpfile << std::scientific << dnam.get_energy() << "\t" << 0. << "\t" << dnam.get_avg_dist() << endl;
	
	ofstream logfile("logfile.txt");
	clock_t starttime = clock();
	for(unsigned k = 0; k < 1000; ++k)
	{
	  logfile << "Started round " << k << endl;
	  	double accepted = dnam.run(1000000);
		  logfile << "Finished run " << k << endl;
// 		dnam.run(10);
//		cout << std::scientific << dnam.get_energy() << "\t" << dnam.get_acceptance_ratio() << "\t" << dnam.get_avg_dist() << endl;
		tmpfile << std::scientific << dnam.get_energy() << "\t" << accepted << "\t" << dnam.get_avg_dist() << endl;
		logfile << "Post first output " << k << endl;
		
		string rname = string("/imports/bott.work/loeb/dna_metropolis/res_res_") + boost::lexical_cast<string>(k) + string(".txt");
		ofstream rofile(rname.c_str());
		logfile << "Second output file opened " << k << endl;
		rofile << dnam;
		rofile.close();
		logfile << "Post second output " << k << endl;
		
		string dname = string("/imports/bott.work/loeb/dna_metropolis/res_distances_") + boost::lexical_cast<string>(k) + string(".txt");
		ofstream distfile(dname.c_str());
		logfile << "Third output file opened " << k << endl;
		for(unsigned i = 0; i < 5000; ++i)
		{
		  logfile << "\t" << i;
			double rv = dnam.inter_point_distance(0., 1.e4*i);
			distfile << 1.e4*i << "\t" << rv << endl;
		}
		distfile.close();
		logfile << endl << "Post third output " << endl << endl;
	}
	clock_t endtime = clock();
	
// 	vector<pair<double, double> > tres = dnam.get_distance_pairs(1000, 3.e6, max_position);
	
	
	logfile << "Execution time: " << (endtime-starttime)/CLOCKS_PER_SEC << endl;
	logfile.close();
	tmpfile.close();
// 	ofstream outfile("results.txt");
// 	outfile << dnam;
// 	outfile.close();
	return 1;
}
