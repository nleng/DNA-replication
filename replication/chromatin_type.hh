
#ifndef CHROMATIN_TYPE_HH
#define CHROMATIN_TYPE_HH

#include <vector>
#include <utility>

using std::vector;
using std::pair;

struct chromatin_type
{
	chromatin_type()
	: fork_speed(28.), base(.1), sigma(1.e6), induced_limit(.1)
	{}
	
	double fork_speed;
	double base;
	double sigma;
	double induced_limit;
	
	chromatin_type& operator= (const chromatin_type &right)
	{
		fork_speed = right.fork_speed;
		base = right.base;
		sigma = right.sigma;
		induced_limit = right.induced_limit;
		return *this;
	}
};

#endif // CHROMATIN_TYPE_HH
