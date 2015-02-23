
#ifndef CHANGE_PARAMETER_HH
#define CHANGE_PARAMETER_HH

enum cps {CHANGE_REPLI_SPEED, CHANGE_BASE, CHANGE_SIGMA};

struct change_parameter
{
	long double time;
	long double value[2];
	size_t chromatin;
	cps command;
	
	change_parameter()
	: time(0.), chromatin(0), command(CHANGE_BASE)
	{}
	
	change_parameter(const long double inval)
	: time(inval), chromatin(0), command(CHANGE_BASE)
	{}
	
	bool operator< (const change_parameter &right) const
	{
		return time < right.time;
	}
	
	change_parameter& operator= (const change_parameter &right)
	{
		time = right.time;
		value[0] = right.value[0];
		value[1] = right.value[1];
		chromatin = right.chromatin;
		command = right.command;
		return *this;
	}
};


#endif // CHANGE_PARAMETER_HH
