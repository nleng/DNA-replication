// simtools - A C++ toolset for stochastic network dynamics models.
// 
// Copyright © 2010-2012 Daniel Löb <daniel@zombiepiratesfromspace.eu>
// 
// simtools is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// simtools is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with simtools.  If not, see <http://www.gnu.org/licenses/>.


// The binary_sum_search tree is a special kind of search tree in that it
// doesn't only contain entries but also at each node holds sums of all
// entries behind that node.
// This makes it very convenient for usage cases such as the Gillespie
// algorithm with many different kinds of events.

#ifndef BINARY_SUM_SEARCH_HH
#define BINARY_SUM_SEARCH_HH

#include <vector>
#include <ostream>
#include "twoe.hh"

using std::vector;
using std::ostream;
using std::endl;
using std::scientific;


namespace simtools
{
	template<typename data_t> class binary_sum_search
	{
		public:
			binary_sum_search(const vector<data_t> &input);
			~binary_sum_search();
			
			data_t get_propensity() const // Returns the probability value at the root of the tree.
			{
				return data[0];
			}
			
			unsigned long get_index(const data_t &randnum) const;
			
			template<typename d_t> friend ostream& operator<<(ostream &the_stream, const binary_sum_search<d_t> &bis);
			
			void rebuild(const unsigned long &position, const data_t &value);	// Change the value in leaf Nr. position to value value and
											// perform all the neccessary modifications to the tree.
			long unsigned get_depth() const
			{
				return depth;
			}
			
			data_t get_value_from_index(const long unsigned &position) const;
			
		protected:
			void full_rebuild();	// Rebuild the whole tree, assuming that only the values in the leaves exist.
			long unsigned depth;	// How many layers deep the tree is. Starts counting at one.
			data_t *data;
	};


	// Retrieve the index of that entry corresponding to the random number randnum, for which
	// 0 <= randnum <= data[0] must hold.
	template<typename data_t> unsigned long binary_sum_search<data_t>::get_index(const data_t &randnum) const
	{
		if ((randnum < 0.) or (randnum > data[0]))
			throw "Invalid random number passed to get_index.";
		
		unsigned long index = 0;
		long unsigned level = 1;
		data_t temprand = randnum;
		
		while (level < depth)
		{
			long unsigned relind = 2*index;
			long unsigned leftind = relind + twoe(level)-1;
			
			if(temprand > data[leftind])
			{
				index = relind+1;
				temprand -= data[leftind];
			}
			else
			{
				index = relind;
			}
			level++;
		}
		return index;
	}

	template<typename data_t> binary_sum_search<data_t>::binary_sum_search(const vector<data_t> &input)
	{
		// Up front, determine the depth of this search tree.
		long unsigned depth_counter = 1;
		unsigned long n_entries = 1;
		while (n_entries < input.size())
		{
			depth_counter++;
			n_entries *= 2;
		}
		depth = depth_counter;
		
		// Allocate data block and fill tree leaves.
		data = new data_t[n_entries*2-1];
		
		data_t *current_pos = data + n_entries-1;
		for(typename vector<data_t>::const_iterator p = input.begin(); p != input.end(); ++p, ++current_pos)
		{
			*current_pos = *p;
		}
		
		// Fill all remaining leaves with zeros.
		while (current_pos < data + n_entries*2-1)
		{
			*current_pos = 0.;
			++current_pos;
		}
		
		full_rebuild();
	}

	template<typename data_t> data_t binary_sum_search<data_t>::get_value_from_index(const long unsigned &position) const
	{
		long unsigned offset = twoe(depth-1)-1;
		return data[offset+position];
	}

	template<typename data_t> void binary_sum_search<data_t>::rebuild(const unsigned long &position, const data_t &value)
	{
		long unsigned temp_depth = depth;
		unsigned long tempos = position;
		long unsigned offset = twoe(depth-1)-1;
		data[offset+position] = value;
		
		while (temp_depth > 1)
		{
			long unsigned oldof = offset;
			temp_depth--;
			offset = (offset+1)/2 -1;
			tempos = tempos/2;
			data[offset+tempos] = data[oldof + 2*tempos] + data[oldof + 2*tempos+1];
		}
	}

	template<typename data_t> void binary_sum_search<data_t>::full_rebuild()
	{
		long unsigned temp_depth = depth;
		while (temp_depth > 1)
		{
			long unsigned offset = twoe(temp_depth-1)-1;
			long unsigned ominusone = twoe(temp_depth-2)-1;
			long unsigned n_entries = twoe(temp_depth-2);
			
			for(long unsigned i = 0; i < n_entries; i++)
			{
				data[ominusone + i] = data[offset + 2*i] + data[offset + 2*i+1];
			}
			
			temp_depth--;
		}
	}

	template<typename data_t> binary_sum_search<data_t>::~binary_sum_search()
	{
		delete[] data;
	}


	template <typename d_t> ostream& operator<<(ostream& the_stream, const binary_sum_search<d_t> &bis)
	{
		long unsigned depth = bis.depth;
		
		// Build up the list of vertical inter-word distances.
		
		vector<long unsigned> idists;
		
		for (long unsigned i = depth; i > 0; --i)
		{
			long unsigned interdist = twoe(i-1);
			idists.push_back(interdist);
		}
		
		// Now step through the data and in each line check, which branches are to be written.
		
		long unsigned maxv = twoe(depth-1);
		for(long unsigned i = 0; i < maxv; ++i)
		{
			for (long unsigned j = 0; j < depth; ++j)
			{
				if (i%idists[j] == 0)
				{
					long unsigned tmp = depth;
					long unsigned t_i = i;
					while (tmp > j+1)
					{
						t_i = t_i/2;
						tmp--;
					}
					long unsigned offset = twoe(j)-1;
					
					the_stream << scientific << bis.data[offset+t_i] << "\t";
				}
				else
				{
					the_stream << "\t" << "\t";
				}
			}
			the_stream << endl;
		}
		return the_stream;
	}
} // End of namespace simtools.

#endif // BINARY_SUM_SEARCH_HH