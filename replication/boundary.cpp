
#include "boundary.hh"

boundary::boundary(double inpos, bool icb, size_t inleft, size_t inright, int ilc, int irc)
	: position(inpos), left(inleft), right(inright), left_chromosome(ilc), right_chromosome(irc), is_chromosome_boundary(icb)
{
}

boundary::boundary(const boundary &rightb)
{
	position = rightb.position;
	left = rightb.left;
	right = rightb.right;
	left_chromosome = rightb.left_chromosome;
	right_chromosome = rightb.right_chromosome;
	is_chromosome_boundary = rightb.is_chromosome_boundary;
}

bool boundary::operator< (const boundary &right) const
{
	return position < right.position;
}


bool boundary::operator== (const boundary &right) const
{
	return position == right.position;
}
