
#include "fork.hh"

fork_t::fork_t(double inpos, double intime, unsigned int inoid, unsigned int inid, unsigned int inchrom, int inchromosome)
	: position(inpos), time(intime), id(inid), anni(0), origin_id(inoid), chromatin(inchrom), start_time(intime), chromosome(inchromosome)
{
}
