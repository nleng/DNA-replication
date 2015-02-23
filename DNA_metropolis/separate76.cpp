
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
#include <boost/filesystem.hpp>

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
using boost::lexical_cast;
// using boost::filesystem::copy_file;
using boost::filesystem::create_directories;

int main(void)
{
	unsigned  calculation_number = 0;
	fstream hn_file("/etc/hostname",std::ios::in);
	char hn[128];
	hn_file.getline(hn,127);
	hn_file.close();
	string outdir = string("./");
	
	string fulldir = outdir + string("dna_metropolis/separate_")+lexical_cast<string>(calculation_number);
	if (not create_directories(fulldir.c_str()))
	{
		cerr << "Failed to access work directory." << endl;
	}
	
	string mydest = fulldir + string("/main.cpp");
// 	copy_file("/home/loeb/phd/programs/dna_metropolis/main_twozones.cpp", mydest ,boost::filesystem::copy_option::overwrite_if_exists);
	
// 	fstream infile("/home/loeb/phd/programs/dna_metropolis/chromosome_patchlists.txt");
	
	vector<list<pair<size_t, double> > > band_data;
	
	for (unsigned ik = 0; ik < 76; ++ik)
	{
		
		list<pair<size_t, double> > t_bands;
		t_bands.push_back(pair<size_t, double>(0,double((1+ik)*1.2e5)));
		band_data.push_back(t_bands);
	}
	
// 	while(infile)
// 	{
//                 list<pair<size_t, double> > t_bands;
//                 size_t chromtype = 0;
//                 string t_line;
//                 getline(infile, t_line);
//                 istringstream line(t_line);
//                 
//                 line.ignore(2,'\n');
//                 line >> chromtype;
//                 
// 		double curpos = 0.;
//                 while(not line.eof())
//                 {
//                         double tmp;
//                         line >> tmp;
// 			curpos += tmp;
// //                      cout << line.good() << "\t" << line.fail() << "\t" << line.eof() << "\t" << line.bad() << endl;
//                         t_bands.push_back(pair<size_t, double>(chromtype,curpos));
// 			chromtype = (chromtype == 1)? 0:1;
//                 }
//                 
// //                 pair<size_t, list<double> > rpair(chromtype, t_bands);
//                 band_data.push_back(t_bands);
//         }

	
	vector<unsigned> n_chromosomes = {76};
	vector<unsigned> n_beads = {1};
	vector<double> real_ends = {0.4};
	
	unsigned beadnr = 0;
	unsigned cntr = 0;
	for(vector<unsigned>::iterator p = n_chromosomes.begin(); p != n_chromosomes.end(); ++p)
	{
		beadnr += (*p)*(n_beads[cntr] + 1); // Plus one because of the skipped bead in between two chromosomes.
		cntr++;
	}
	dna_metropolis dnam(beadnr, 0., 0., 0, 1.e3, 290);
	
	dnam.use_loops = false;
	dnam.use_gravity = true;
// 	dnam.start_big_random = true;
	dnam.init(n_chromosomes, n_beads, real_ends, band_data, 0, 0, 0);
	
	string fname1 = fulldir + string("separate_")+lexical_cast<string>(calculation_number) + string("_tmpfile.txt");
	ofstream tmpfile(fname1.c_str());
// 		cout << std::scientific << dnam.get_energy() << "\t" << dnam.get_acceptance_ratio() << "\t" << dnam.get_avg_dist() << endl;
	tmpfile << std::scientific << dnam.get_energy() << "\t" << 0. << "\t" << dnam.get_avg_dist() << endl;
	string fname2 = fulldir + string("separate")+lexical_cast<string>(calculation_number) + string("_logfile.txt");
	string fname3 = fulldir + string("separate")+lexical_cast<string>(calculation_number) + string("_centers.txt");
	ofstream logfile(fname2.c_str());
	ofstream chromfile(fname3.c_str());
	clock_t starttime = clock();
	for(unsigned k = 0; k < 1000; ++k)
	{
		logfile << "Started round " << k << endl;
		double accepted = dnam.run(10000);
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
		
 		string rname = fulldir + string("/res_res_") + boost::lexical_cast<string>(k) + string(".txt");
		ofstream rofile(rname.c_str());
		//	logfile << "Second output file opened " << k << endl;
		rofile << dnam;
		rofile.close();
		//logfile << "Post second output " << k << endl;
		//		string dname = string("results/ranlo/res_distances_") + boost::lexical_cast<string>(k) + string(".txt");
 		string dname = fulldir + string("/res_distances_") + boost::lexical_cast<string>(k) + string(".txt");
		ofstream distfile(dname.c_str());
		//logfile << "Third output file opened " << k << endl;
		for(unsigned i = 0; i < 20000; ++i)
		{
		  //logfile << "\t" << i;
			double rv = dnam.inter_point_distance(0., 1.e4*i);
			distfile << 5.e4*i << "\t" << rv << endl;
		}
		distfile.close();
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
