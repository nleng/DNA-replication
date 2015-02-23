#ifndef ANNIHILATION_HH
#define ANNIHILATION_HH

#include "fork.hh"
#include "boundary.hh"
#include "simtools/red_black_tree.hh"

using simtools::red_black_tree;

// enum ist aufzaehlungstyp, hier: art der annihilation, fork+fork oder fork+boundary, aber was ist solo removal?
enum anni_type {FORK_FORK, FORK_BOUNDARY, SOLO_REMOVAL};

struct annihilation
{
	annihilation(double inval = 0.)
	: endtime(inval), at(FORK_FORK), first_fork(NULL), second_fork(NULL), other_boundary(NULL), left(true)
	{}
	
	annihilation(const annihilation &input);
	bool operator< (const annihilation &right) const;
	annihilation& operator= (const annihilation &right);
	bool hasequal_fork(const annihilation &right);
	
	long double endtime;
	anni_type at;
	
	red_black_tree<fork_t>::iterator first_fork;
	red_black_tree<fork_t>::iterator second_fork;
	red_black_tree<boundary>::iterator other_boundary;
	bool left;
};

#endif // ANNIHILATION_HH
