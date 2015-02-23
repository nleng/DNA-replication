
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
using boost::filesystem::create_directories;


int main(void)
{
	fstream infile("chromosome_patchlists.txt");
        
        vector<list<pair<size_t, double> > > band_data;
        
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
                while(not line.eof())
                {
                        double tmp;
                        line >> tmp;
			curpos += tmp;
//                      cout << line.good() << "\t" << line.fail() << "\t" << line.eof() << "\t" << line.bad() << endl;
                        t_bands.push_back(pair<size_t, double>(chromtype,curpos));
			chromtype = (chromtype == 1)? 0:1;
                }
                
//                 pair<size_t, list<double> > rpair(chromtype, t_bands);
                band_data.push_back(t_bands);
        }
	
	vector<unsigned> n_chromosomes = {5,3,3,2,6,3,5,3,5,3,3,3,3,3,3,3,4,2,3,3,3,3,2};
	vector<unsigned> n_beads = {2259, 2426, 1987, 1912, 1795, 1707, 1581, 1462, 1169, 1348, 1328, 1320, 1121, 1047, 983, 875, 786, 754, 614, 620, 446, 488, 1540, 271};
	vector<double> real_ends = {225949719.0, 242651149.0, 198701828.0, 191273063.0, 179557866.0, 170799992.0, 158121424.0, 146274826.0, 116973252.0, 134874737.0, 132852384.0, 132049534.0, 112142981.0, 104768586.0, 98338916.0, 87527254.0, 78674742.0, 75417153.0, 61411651.0, 62035964.0, 44644324.0, 48891433.0, 154013754.0, 27100000.0};
// 	vector<double> real_ends = {247249719.0, 242951149.0, 199501827.0, 191273063.0, 180857866.0, 170899992.0, 158821424.0, 146274826.0, 140273252.0, 135374737.0, 134452384.0, 132349534.0, 114142980.0, 106368585.0, 100338915.0, 88827254.0, 78774742.0, 76117153.0, 63811651.0, 62435964.0, 46944323.0, 49691432.0, 154913754.0}

	unsigned calculation_number = 100;
	fstream hn_file("/etc/hostname",std::ios::in);
        char hn[128];
        hn_file.getline(hn,127);
        hn_file.close();
	string outdir = string("/imports/") + string(hn) + string(".work/loeb/");
        string fulldir = outdir + string("dna_metropolis/loops_")+lexical_cast<string>(calculation_number);
	if (not create_directories(fulldir.c_str()))
	{
	    cerr << "Failed to access work directory." << endl;
	}	
	
	unsigned beadnr = 0;
	unsigned cntr = 0;
	for(vector<unsigned>::iterator p = n_chromosomes.begin(); p != n_chromosomes.end(); ++p)
	{
		beadnr += (*p)*n_beads[cntr] + 1; // Plus one because of the skipped bead in between two chromosomes.
		cntr++;
	}
	dna_metropolis dnam(beadnr, 8.e-8, 8.e-8, 5.e-7, 0.0002, 290);
	dnam.ellipsis_y = 5000;	

	unsigned n_connections = 5000;
	dnam.use_loops = true;
	dnam.use_gravity = false;
// 	dnam.start_big_random = true;
	dnam.init(n_chromosomes, n_beads, real_ends, band_data, n_connections*7/11, n_connections*3/11, n_connections/11, false);
	string fname2 = string("loops_")+lexical_cast<string>(calculation_number) + string("_logfile.txt");
	//        string fname3 = string("loops_")+lexical_cast<string>(calculation_number) + string("_centers.txt");
	string fname1 = string("loops_")+lexical_cast<string>(calculation_number) + string("_tmpfile.txt");
        ofstream tmpfile(fname1.c_str());

	tmpfile << std::scientific << dnam.get_energy() << "\t" << 0. << "\t" << dnam.get_avg_dist() << endl;
	
	ofstream logfile(fname2.c_str());
	clock_t starttime = clock();
	for(unsigned k = 0; k < 1000; ++k)
	{
	  logfile << "Started round " << k << endl;
	  	double accepted = dnam.run(10000000);
		//	  logfile << "Finished run " << k << endl;
// 		dnam.run(10);
//		cout << std::scientific << dnam.get_energy() << "\t" << dnam.get_acceptance_ratio() << "\t" << dnam.get_avg_dist() << endl;
		tmpfile << std::scientific << dnam.get_energy() << "\t" << accepted << "\t" << dnam.get_avg_dist() << endl;
		//logfile << "Post first output " << k << endl;

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
			distfile << 1.e4*i << "\t" << rv << endl;
		}
		distfile.close();
		//logfile << endl << "Post third output " << endl << endl;
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
