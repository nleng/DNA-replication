
#include <vector>
#include <iostream>
#include <fstream>
#include "dna_metropolis.hh"
#include "get_plotpoints.hh"
#include <utility>
#include <string>
#include <sstream>
#include <list>
#include <set>
#include <boost/lexical_cast.hpp>

using namespace boost::python;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::pair;
using std::istringstream;
using std::fstream;

list get_3d_points(list in_forks, str beadfilename)
{
	fstream infile("chromosome_patchlists.txt");
        
        vector<std::list<pair<size_t, double> > > band_data;
        
        while(infile)
        {
                std::list<pair<size_t, double> > t_bands;
                size_t chromtype = 0;
                std::string t_line;
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
	vector<unsigned> n_beads = {2260, 2427, 1988, 1913, 1796, 1708, 1582, 1463, 1170, 1349, 1329, 1321, 1122, 1048, 984, 876, 787, 755, 615, 621, 447, 489, 1541, 272};
	vector<double> real_ends = {225949719.0, 242651149.0, 198701828.0, 191273063.0, 179557866.0, 170799992.0, 158121424.0, 146274826.0, 116973252.0, 134874737.0, 132852384.0, 132049534.0, 112142981.0, 104768586.0, 98338916.0, 87527254.0, 78674742.0, 75417153.0, 61411651.0, 62035964.0, 44644324.0, 48891433.0, 154013754.0, 27100000.0};
	
	unsigned beadnr = 0;
	unsigned cntr = 0;
	for(vector<unsigned>::iterator p = n_chromosomes.begin(); p != n_chromosomes.end(); ++p)
	{
		beadnr += (*p)*n_beads[cntr] + 1; // Plus one because of the skipped bead in between two chromosomes.
		cntr++;
	}
// 	beadnr--; // Last one needs no skipped bead. // Cancelled out by additional first bead.
	
	dna_metropolis dnam(beadnr, 6.e-5, 6.e-5, 6.e-5, 4., 290.);
	
	unsigned n_connections = 5000;
	dnam.use_loops = true;
	dnam.use_gravity = false;
// 	dnam.start_big_random = true;
	dnam.init(n_chromosomes, n_beads, real_ends, band_data, n_connections*7/11, n_connections*3/11, n_connections/11, true);
	
	std::string bfname = extract<std::string>(beadfilename);
	fstream coordfile(bfname.c_str());
	
        vector<double> coordinates;
	coordinates.reserve(120000);
	
        while(coordfile)
        {
		std::string t_line;
		getline(coordfile, t_line);
		istringstream line(t_line);
		
		while(not line.eof())
                {
			double tmp;
			line >> tmp;
			coordinates.push_back(tmp);
                }
        }
        coordfile.close();
        
	dnam.load_coordinates(coordinates);
	
	std::set<double> forks;
	
	for(unsigned int i = 0; i < len(in_forks); ++i)
	{
		forks.insert(extract<double>(in_forks[i]));
		
	}
	
// 	unsigned flen = forks.size();
	
	list eu_results;
	list het_results;
	for(std::set<double>::iterator p = forks.begin(); p != forks.end(); ++p)
	{
		vertex tv = dnam.dna_to_position(*p);
		size_t chromt = dnam.dna_to_chromatin_type(*p);
		list tlist;
		tlist.append(object(tv.x));
		tlist.append(object(tv.y));
		tlist.append(object(tv.z));
		
		if (chromt)
			het_results.append(tlist);
		else
			eu_results.append(tlist);
	}
	
	list results;
	results.append(eu_results);
	results.append(het_results);
	return results;
}

// Expose all classes and functions to python

BOOST_PYTHON_MODULE(libreplihelpers)
{
	def("get_3d_points", get_3d_points);
}