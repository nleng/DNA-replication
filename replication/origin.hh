
#ifndef ORIGIN_HH
#define ORIGIN_HH

struct origin
{
	origin(double ipos, size_t ichro, int inchromosome)
	: position(ipos), chromosome(inchromosome), chromatin(ichro)
	{}
	
	bool operator< (const origin &right) const
	{
		return position < right.position;
	}
	
	bool operator== (const origin &right) const
	{
		// Only take the position into account here!
		if (position == right.position)
			return true;
		return false;
	}
	
	double position;
	int chromosome;
	size_t chromatin;	// Identifies the chromatin type that this origin is in.
};

#endif // ORIGIN_HH
