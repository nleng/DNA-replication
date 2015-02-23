#include "twoe.hh"

unsigned long twoe(const size_t &in)
{
	if (in > 63)
		throw "twoe: Cannot store number greater than 2exp64 -1 in unsigned long.";
	
	unsigned long ret = 1;
	return (ret << in);
}
