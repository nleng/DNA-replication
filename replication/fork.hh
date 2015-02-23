
#ifndef FORK_HH
#define FORK_HH

struct fork_t
{
	fork_t(double inpos, double intime, unsigned int inoid, unsigned int inid, unsigned int inchrom, int inchromosome);
	
	bool operator< (const fork_t &right) const
	{
		return position < right.position;
	}
	
	double position;
	double time;
	unsigned int id;	// Unambiguous ID of this fork.
	unsigned int anni;	// Id of the annihilation event.
	unsigned int origin_id;	// Id of the original firing event that fork came from.
	unsigned int chromatin;	// Chromatin type that this fork is currently in.
	double start_time;
	int chromosome;
};

#endif // FORK_HH
