

#include <vector>
#include <set>
#include <list>
#include <map>
#include <utility>
#include <tr1/random>
#include <boost/generator_iterator.hpp>
// #include <ostream>
#include <iostream>
#include <stdlib.h>     /* srand, rand */

struct vertex
{
	double x;
	double y;
	double z;
};

struct chromosome
{
	chromosome()
	: end_pos(0.), start_pos(0.), index(0)
	{}
	
	bool operator< (const chromosome &right) const
	{
		return end_pos < right.end_pos;
	}
	
	unsigned end_pos;
	unsigned start_pos;
	unsigned index;
};

class dna_metropolis
{
	public:
		dna_metropolis(unsigned int n_beads, double i_kappa_het, double i_kappa_eu, double i_kappa_fac, double i_kappa_2, double i_kappa_3, double i_kappa_fac_rep, double temperature);
		double run(unsigned n_steps);
		double get_energy() const {return U;}
// 		double get_acceptance_ratio() const {return double(n_accepted)/double(n_cur_step);}
		double get_avg_dist() const;
		double inter_point_distance(double posa, double posb);
		vertex dna_to_position(double dna_pos) const;
		size_t dna_to_chromatin_type(double dna_pos) const;
		std::vector<std::pair<double, double> > get_distance_pairs(unsigned n_pairs, double max_dist, double max_position);
		void prepare_chromosomes();
		void init(const std::vector<unsigned> &n_chromosomes, const std::vector<unsigned> &chromo_beads, const std::vector<double> &real_ends, const std::vector<std::list<std::pair<size_t, double> > > &band_data, unsigned eu_cons, unsigned het_cons, unsigned fac_cons, unsigned inter_cons, bool truly_random_connections);
		void load_coordinates(const std::vector<double> &in_coords);
		
		friend std::ostream& operator<<(std::ostream &the_stream, const dna_metropolis &dnam);
		
		double ellipsis_x;
		double ellipsis_y;
		double ellipsis_z;
	protected:
		void build_chromatin();
		bool isinside(const double *p) const;
		double en_calc(const double *left, const double *right, const unsigned &binf) const; 	// Helper function.
		double fac_repuls(const double*, const unsigned int&) const;
// 		double en_calc(const double *left, const double *middle, const double *right) const;
		double inter_loop_en_calc(const double *left, const double *right) const;
// 		double inter_loop_en_calc(const double *left, const double *middle, const double *right) const;
		double center_delta_en_calc(const double *coords, unsigned index, unsigned chromo_index) const;
		double calc_centers_energy() const;
		double calc_energy();				// Initial total energy calculation.
		bool loop_relevant(unsigned index) const {return (bead_info[index]>>9)%2;};
		
		typedef std::tr1::mt19937 engine_t;
		typedef std::tr1::uniform_int<size_t> uniform_int_t;
		typedef std::tr1::uniform_real<long double> uniform_real_t;
		typedef std::tr1::variate_generator<engine_t&, uniform_int_t> uni_int_gen;
		typedef std::tr1::variate_generator<engine_t&, uniform_real_t > vg_real_t;
		typedef boost::generator_iterator<vg_real_t > real_it_t;
		typedef boost::generator_iterator<uni_int_gen> int_it_t;
		// Spring constant parameters
		
		unsigned n_beads;
		double kappa_het;
		double kappa_eu;
		double kappa_fac;
		double kappa_2;
		double kappa_3;
		double kappa_fac_rep;
		
		double beta;
		
		constexpr static double boltzmann_constant = 1.380648e-23;
		// The unit of measurement here is nanometers.
		constexpr static double init_dist = 200.;	// Inter-bead distance used for initial placement of beads.
		constexpr static double delta_jump = 10.;
		constexpr static double dna_unit = 1.e5;
		constexpr static double pi = 3.14159265;
		constexpr static unsigned connection_indicator = 131071; // If this value is set for a connection, then there is no connection.
		
		double U;	// Current system energy value.
		unsigned long n_cur_step;
		unsigned long n_accepted;
		
		std::vector<double> coordinates;
		std::vector<unsigned int> bead_info;	// Contains one entry for every bead indicating if it starts/ends a chromosome and if it is a loop connection.
		std::set<chromosome> chromosomes;
		std::vector<double> chromo_positions;
		std::vector<double> chromo_weights;
	public:
// 		std::set<unsigned int> chromosome_starters;
// 		std::set<unsigned int> chromosome_enders;
		std::map<double, std::pair<unsigned,double> > position_offsets;
		bool use_loops;
		bool use_gravity;
  		bool start_big_random;
		bool nonucleus;
		bool with_nucleoli;
		double surface_rad;
		double hetFactor;
		double nucleolusRad;
		double nucleolusRad2;
		const std::vector<double> nucleolusPos;
		const std::vector<double> nucleolusPos2;
		unsigned long get_step_number() const { return n_cur_step;}
		std::vector<double> get_chromosome_centers() const {return chromo_positions;};
	protected:	
		engine_t rng_engine;
		uniform_real_t rng_real;
		vg_real_t vg_real;
		real_it_t real_it;
		
		uniform_int_t rng_int;
		uni_int_gen vg_int;
		int_it_t int_it;
};

std::ostream& operator<<(std::ostream &the_stream, const dna_metropolis &dnam);
