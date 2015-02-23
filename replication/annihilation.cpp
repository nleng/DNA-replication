
#include "fork.hh"
#include "annihilation.hh"

// annihilation::annihilation (double inval = 0.)
// 	: endtime(inval), at(FORK_FORK), first_fork(NULL), second_fork(NULL), other_boundary(NULL)
// {}

annihilation::annihilation(const annihilation &input)
	: endtime(input.endtime), at(input.at), first_fork(input.first_fork), second_fork(input.second_fork), other_boundary(input.other_boundary), left(input.left)
{
}

bool annihilation::operator< (const annihilation &right) const
{
	return endtime < right.endtime;
}

bool annihilation::hasequal_fork(const annihilation &right)
{
	if ((first_fork == right.first_fork) or (second_fork == right.first_fork) or (first_fork == right.second_fork) or ((second_fork == right.second_fork) and (second_fork != NULL)))
		return true;
	return false;
}

annihilation& annihilation::operator= (const annihilation &right)
{
	endtime = right.endtime;
	at = right.at;
	
	first_fork = right.first_fork;
	second_fork = right.second_fork;
	other_boundary = right.other_boundary;
	left = right.left;
	
	return *this;
}
