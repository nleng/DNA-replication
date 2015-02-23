

#ifndef BOUNDARY_HH
#define BOUNDARY_HH

#include <stddef.h>

struct boundary
{
	boundary(double inpos, bool icb, size_t inleft, size_t inright, int ilc, int irc);
	boundary(const boundary &rightb);
	bool operator< (const boundary &right) const;
	bool operator== (const boundary &right) const;
	double position;
	
	// Left and right chromatin identifiers.
	size_t left;
	size_t right;
	
	int left_chromosome;
	int right_chromosome;
	
	// Is this a hard (i.e. chromosome) boundary?
	bool is_chromosome_boundary;
};

#endif // BOUNDARY_HH
