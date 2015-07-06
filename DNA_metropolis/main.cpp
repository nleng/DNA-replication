
#include <vector>
#include <iostream>
#include <fstream>
#include "dna_metropolis.hh"
#include <utility>
#include <ctime>
#include <string>
#include <sstream>
#include <list>
#include <boost/lexical_cast.hpp>

using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::pair;
using std::clock;
using std::string;
using std::istringstream;
using std::list;
using std::fstream;

int main(void)
{
    fstream infile2("chromatinType.txt");
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

	fstream infile("chromosome_patchlists_inactiveX.txt");
	vector<list<pair<size_t, double> > > band_data;
	int chromCounter = 0;
	while(infile)
	{
		list<pair<size_t, double> > t_bands;
		size_t chromtype = 0;
		string t_line;
		getline(infile, t_line);
		istringstream line(t_line);

		line.ignore(2,'\n');
		line >> chromtype;

		double curpos = 0.;
		int zoneCounter = 0;
		while(not line.eof())
		{
			double tmp;
			line >> tmp;
			curpos += tmp;
			chromtype = chromType_data[chromCounter][zoneCounter];
			t_bands.push_back(pair<size_t, double>(chromtype,curpos));
			++zoneCounter;
		}
		band_data.push_back(t_bands);
		++chromCounter;
	}
	
	vector<unsigned> n_chromosomes = {5,3,3,2,6,3,5,3,5,3,3,3,3,3,3,3,4,2,3,3,3,3,1,1};
	vector<unsigned> n_beads = {2473, 2430, 1996, 1913, 1809, 1709, 1589, 1463, 1403, 1354, 1345, 1324, 1142, 1064, 1004, 889, 788, 762, 639, 625, 470, 497, 1550, 1550};
	vector<double> real_ends = {247249719.0, 242951149.0, 199501827.0, 191273063.0, 180857866.0, 170899992.0, 158821424.0, 146274826.0, 140273252.0, 135374737.0, 134452384.0, 132349534.0, 114142980.0, 106368585.0, 100338915.0, 88827254.0, 78774742.0, 76117153.0, 63811651.0, 62435964.0, 46944323.0, 49691432.0, 154913754.0, 154913754.0};
	
	unsigned beadnr = 0;
	unsigned cntr = 0;
	for(vector<unsigned>::iterator p = n_chromosomes.begin(); p != n_chromosomes.end(); ++p)
	{
		beadnr += (*p)*(n_beads[cntr] + 1); // Plus one because of the skipped bead in between two chromosomes.
		cntr++;
	}
	// kappa_2: for the loop constant
	// kappa_3: for the repulsion of the chromatin-center
	// n_beads, kappa_het, kappa_eu, kappa_fac, kappa_2, kappa_3, kappa_fac_rep, temperature)
	string outdir = ("/my/path/");
	dna_metropolis dnam(beadnr, 3.e-6, 1.e-8, 5.e-7, 5.e-7, 1.e-4,  30., 290);
	dnam.ellipsis_y = 5000;
	
	unsigned n_connections = 5000;
	dnam.use_loops = true;
	dnam.use_gravity = true;
// 	dnam.start_big_random = true;
	// number of connections intra and inta chromatin type
	dnam.init(n_chromosomes, n_beads, real_ends, band_data, n_connections/5, n_connections/5, n_connections/5, n_connections/5*2, true);
// 	dnam.init(n_chromosomes, n_beads, real_ends, band_data, n_connections*3/16, n_connections*7/16, n_connections*5/16, n_connections/16, false);

// 	string outdir = ("/imports/debora.work/nicor/dnaMetropolisResults/x0/");
// 	string outdir = ("/home/corni/aaReplication/dnaSim/outFile/");
	ofstream tmpfile((outdir + string("tmpfile.txt")).c_str());
// 		cout << std::scientific << dnam.get_energy() << "\t" << dnam.get_acceptance_ratio() << "\t" << dnam.get_avg_dist() << endl;
	tmpfile << std::scientific << dnam.get_energy() << "\t" << 0. << "\t" << dnam.get_avg_dist() << endl;
	
	ofstream logfile((outdir + string("logfile.txt")).c_str());
	ofstream chromfile((outdir + string("centers.txt")).c_str());
	clock_t starttime = clock();
	for(unsigned k = 1; k < 10001; ++k)
	{	
		cout  <<  "simNr: "<< k <<  endl;
		logfile << "Started round " << k << endl;
	  	double accepted = dnam.run(10000000);
		//	  logfile << "Finished run " << k << endl;
// 		dnam.run(10);
//		cout << std::scientific << dnam.get_energy() << "\t" << dnam.get_acceptance_ratio() << "\t" << dnam.get_avg_dist() << endl;
		tmpfile << std::scientific << dnam.get_energy() << "\t" << accepted << "\t" << dnam.get_avg_dist() << endl;
		//logfile << "Post first output " << k << endl;
		
		//		string rname = string("results/ranlo/res_res_") + boost::lexical_cast<string>(k) + string(".txt");
		vector<double> ccenters = dnam.get_chromosome_centers();
		for(vector<double>::iterator p = ccenters.begin(); p != ccenters.end(); ++p)
		{
			chromfile << *p << "\t";
		}
		chromfile << endl;
		if(k % 1000 == 0)	// (k % 2000 == 0)
		{
			string rname = outdir + string("res_res_") + boost::lexical_cast<string>(k) + string(".txt");
			ofstream rofile(rname.c_str());
			//	logfile << "Second output file opened " << k << endl;
			rofile << dnam;
			rofile.close();
		}
		//logfile << "Post second output " << k << endl;
		//		string dname = string("results/ranlo/res_distances_") + boost::lexical_cast<string>(k) + string(".txt");
		// 		string dname = string("/imports/hus.work/loeb/dna_metropolis/fullcalc16/res_distances_") + boost::lexical_cast<string>(k) + string(".txt");
		//		ofstream distfile(dname.c_str());
		//logfile << "Third output file opened " << k << endl;
		//		for(unsigned i = 0; i < 20000; ++i)
		//		{
		  //logfile << "\t" << i;
		//			double rv = dnam.inter_point_distance(0., 1.e4*i);
		//			distfile << 1.e4*i << "\t" << rv << endl;
		//		}
		//		distfile.close();
		//logfile << endl << "Post third output " << endl << endl;
	}
	clock_t endtime = clock();
	chromfile.close();
// 	vector<pair<double, double> > tres = dnam.get_distance_pairs(1000, 3.e6, max_position);
	
	
	logfile << "Execution time: " << (endtime-starttime)/CLOCKS_PER_SEC << endl;
	logfile.close();
	tmpfile.close();
// 	ofstream outfile("results.txt");
// 	outfile << dnam;
// 	outfile.close();
	return 1;
}
