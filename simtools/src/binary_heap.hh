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


#ifndef BINARY_HEAP_HH
#define BINARY_HEAP_HH

#include <vector>
#include <map>
#include <ostream>
#include <utility>
#include <iostream>
#include <limits>
#include "twoe.hh"
#include "exceptions.hh"

using std::vector;
using std::ostream;
using std::endl;
using std::scientific;
using std::pair;
using std::map;
using std::cerr;
using std::cout;
using std::endl;
using std::numeric_limits;


namespace simtools
{
	// This is a textbook style implementation of a binary heap. 
	// A typical usage scenario is the storage of time-ordered events.
	template<typename data_t> class binary_heap
	{
		public:
			binary_heap();
			binary_heap(const binary_heap<data_t> &inheap);
			~binary_heap();
			
			pair<data_t, long unsigned> get_smallest() const // Returns the smallest value (located at the heap base)
			{
				if (n_elements == 0)
				{
					pair<data_t, long unsigned> temp(numeric_limits<data_t>::max(),0);
					return temp;
				}
				return data[0];
			}
			
			pair<data_t, long unsigned> remove_smallest();	// Removes the smallest value (located at the heap base) and returns it.
			void remove(const long unsigned &index);
			
			data_t get_value(const long unsigned &index) const;
			
			template<typename d_t> friend ostream& operator<<(ostream &the_stream, const binary_heap<d_t> &bh);
			
			void insert_value(const data_t &input, const long unsigned &index);	// Insert a new value into the tree.
			
			long unsigned get_n_elements() const
			{
				return n_elements;
			}
			
			long unsigned get_depth() const
			{
				return depth;
			}
			
			bool exists(const long unsigned &position) const;
			bool test_consistency() const;
			
			void clear();
		protected:
			long unsigned get_pos_depth(const long unsigned &position) const;
			
			long unsigned depth;	// How many layers deep the tree is. Starts counting at one.
			long unsigned n_elements; // How many date elements are actually filled.
			pair<data_t, long unsigned> *data;
			map<long unsigned, long unsigned> positions;	// Store index and position in the data set for each data point.
	};
	
	template<typename data_t> binary_heap<data_t>::binary_heap(const binary_heap<data_t> &inheap)
	{
		depth = inheap.depth;
		n_elements = inheap.n_elements;
		long unsigned max_size = twoe(depth + 1)-1;
		data = new pair<data_t, long unsigned>[max_size];
		for (unsigned int i = 0; i < n_elements; i++)
		{
			data[i] = inheap.data[i];
		}
		positions = inheap.positions;
	}
	
	template<typename data_t> bool binary_heap<data_t>::exists(const long unsigned &position) const
	{
		typename map<long unsigned, long unsigned>::const_iterator it = positions.find(position);
		if (it == positions.end())
			return false;
		return true;
	}
	
	template<typename data_t> pair<data_t, long unsigned> binary_heap<data_t>::remove_smallest()
	{
		pair<data_t, long unsigned> smallest_value = data[0];
		positions.erase(data[0].second);
		n_elements--;
		
		if (n_elements == 0)
		{
			data[n_elements] = pair<data_t, long unsigned>(data_t(0.),0);
			return smallest_value; // We have to leave the method here, because the next commands will increase the number of positions to at least one.
		}
		
		// Move the last value (formerly of index n_elements-1, now n_elements) to the first position.
		data[0] = data[n_elements];
		data[n_elements] = pair<data_t, long unsigned>(data_t(0.),0);
		positions[data[0].second] = 0;
		
		// Now resolve all the ordering violations this causes.
		long unsigned tempdepth = 1;
		long unsigned curindex = 0;
		long unsigned offset = twoe(tempdepth-1)-1;
		long unsigned factual_depth = (n_elements>0?get_pos_depth(n_elements-1):0);
		
		while (tempdepth < factual_depth)
		{
			long unsigned nextoffset = twoe(tempdepth)-1;
			long unsigned ind1 = curindex*2;
			long unsigned ind2 = ind1 + 1;
			
			if (tempdepth == factual_depth-1)
			{
				if (nextoffset+ind1 >= n_elements)
					break;
			}
			
			if ((data[offset + curindex] > data[nextoffset + ind1]) or ((data[offset + curindex] > data[nextoffset + ind2]) and (nextoffset + ind2 < n_elements)))
			{
				long unsigned oldp = offset + curindex;
				long unsigned newp = nextoffset + ind1;
				// Switch places with the smaller one.
				if ((data[nextoffset + ind1] > data[nextoffset + ind2]) and (nextoffset + ind2 < n_elements))
				{
					newp = nextoffset + ind2;
				}
				
				// Perform the switch.
				positions[data[oldp].second] = newp;
				positions[data[newp].second] = oldp;
				
				pair<data_t, long unsigned> temp = data[newp];
				data[newp] = data[oldp];
				data[oldp] = temp;
				curindex = newp-nextoffset;
			}
			else
			{
				break;
			}
			
			tempdepth++;
			offset = nextoffset;
		}
		
		return smallest_value;
	}

	template<typename data_t> void binary_heap<data_t>::remove(const long unsigned &index)
	{
		if (not (exists(index)))
		{
			throw RuntimeError("Tried to remove non-existing element from binary heap");
		}
		
		// When removing any value from the tree, take care of the ensuing packing violation by inserting the last value and then resolving the ordering violation by moving up and down the tree.
		
		long unsigned dpos = positions[index];
		
		positions.erase(index);
		n_elements--;
		
		if (n_elements == 0)
		{
			data[n_elements] = pair<data_t, long unsigned>(data_t(0.),0);
			return; // We have to leave the method here, because the next commands will increase the number of positions to at least one.
		}
		
		if (dpos == n_elements)
		{
			data[n_elements] = pair<data_t, long unsigned>(data_t(0.),0);
			return;
		}
		
		// Move the last value (formerly of index n_elements-1, now n_elements) to the position of the removed one.
		data[dpos] = data[n_elements];
		data[n_elements] = pair<data_t, long unsigned>(data_t(0.),0);
		positions[data[dpos].second] = dpos;
		
		// Now resolve all the ordering violations this causes. First propagate the inserted value outward, as if this was the deletion of the root value.
		long unsigned tempdepth = get_pos_depth(dpos);
		long unsigned offset = twoe(tempdepth-1)-1;
		long unsigned curindex = dpos - offset;
		long unsigned factual_depth = (n_elements>0?get_pos_depth(n_elements-1):0);
		
		
		while (tempdepth < factual_depth)
		{
			long unsigned nextoffset = twoe(tempdepth)-1;
			long unsigned ind1 = curindex*2;
			long unsigned ind2 = ind1 + 1;
			
			if (tempdepth == factual_depth-1)
			{
				if (nextoffset+ind1 >= n_elements)
					break;
			}
			
			if ((data[offset + curindex] > data[nextoffset + ind1]) or ((data[offset + curindex] > data[nextoffset + ind2]) and (nextoffset + ind2 < n_elements)))
			{
				long unsigned oldp = offset + curindex;
				long unsigned newp = nextoffset + ind1;
				// Switch places with the smaller one.
				if ((data[nextoffset + ind1] > data[nextoffset + ind2]) and (nextoffset + ind2 < n_elements))
				{
					newp = nextoffset + ind2;
				}
				
				// Perform the switch.
				positions[data[oldp].second] = newp;
				positions[data[newp].second] = oldp;
				
				pair<data_t, long unsigned> temp = data[newp];
				data[newp] = data[oldp];
				data[oldp] = temp;
				curindex = newp-nextoffset;
			}
			else
			{
				break;
			}
			
			tempdepth++;
			offset = nextoffset;
		}
		
		// Now take care of the possibility, that the inserted former last value is smaller than the value above its new position. This part mimics a standard insertion procedure.
		
		long unsigned nextof = ((offset+1)/2 > 0 ? (offset+1)/2 : 1) -1;
	// 	unsigned long tempos = curindex;
		unsigned long nextpos = curindex/2;
		data_t curval = data[offset + curindex].first;
		data_t nextval = data[nextof + nextpos].first;
		
		while ((tempdepth > 1) and (curval < nextval))
		{
			long unsigned oldp = offset + curindex;
			long unsigned newp = nextof + nextpos;
			
			// Perform the switch.
			positions[data[oldp].second] = newp;
			positions[data[newp].second] = oldp;
			
			pair<data_t, long unsigned> temp = data[newp];
			data[newp] = data[oldp];
			data[oldp] = temp;
			
			// Update for the next step
			offset = nextof;
			curindex = nextpos;
			nextof = ((nextof+1)/2 > 0 ? (nextof+1)/2 : 1) -1;
			nextpos = nextpos/2;
			nextval = data[nextof + nextpos > 0 ? nextof + nextpos: 0].first;
			tempdepth--;
		}
	}

	template<typename data_t> data_t binary_heap<data_t>::get_value(const long unsigned &index) const
	{
		map<long unsigned, long unsigned>::const_iterator temp = positions.find(index);
		if (temp != positions.end())
			return data[temp->second].first;
		return data_t(-1.);	// This is kind of bad style, in Band communication and all, but so far I have no better way of doing this here.
	}
	template<typename data_t> void binary_heap<data_t>::insert_value(const data_t &input, const long unsigned &index)	// Insert a new value into the tree.
	{
		long unsigned max_n_elems = twoe(depth)-1;
		
		if (n_elements == max_n_elems)
		{
			long unsigned new_max_size = twoe(depth + 1)-1;
			pair<data_t, long unsigned> *temp = new pair<data_t, long unsigned>[new_max_size];
			
			pair<data_t, long unsigned> *p = data;
			for (long unsigned counter = 0; counter < max_n_elems; counter++, ++p)
			{
				temp[counter] = *p;
			}
			delete[] data;
			data = temp;
			depth++;
		}
		
		data[n_elements] = pair<data_t, long unsigned>(input, index);
		positions[index] = n_elements;
		n_elements++;
		
		// Now sort out the possible heap ordering violation.
		long unsigned curdepth = get_pos_depth(n_elements-1);
		
		// So now curdepth holds the depth level which the element is at. Crawl up the heap and perform switches as long as there are bigger values above the current position.
		
		long unsigned offset = twoe(curdepth-1)-1;
		long unsigned nextof = ((offset+1)/2 > 0 ? (offset+1)/2 : 1) -1;
		unsigned long tempos = n_elements -1 -offset;
		unsigned long nextpos = tempos/2;
		data_t curval = data[n_elements-1].first;
		data_t nextval = data[nextof + nextpos].first;
		
		while ((curdepth > 1) and (curval < nextval))
		{
			long unsigned oldp = offset + tempos;
			long unsigned newp = nextof + nextpos;
			
			// Perform the switch.
			positions[data[oldp].second] = newp;
			positions[data[newp].second] = oldp;
			
			pair<data_t, long unsigned> temp = data[newp];
			data[newp] = data[oldp];
			data[oldp] = temp;
			
			// Update for the next step
			offset = nextof;
			tempos = nextpos;
			nextof = ((nextof+1)/2 > 0 ? (nextof+1)/2 : 1) -1;
			nextpos = nextpos/2;
			nextval = data[nextof + nextpos > 0 ? nextof + nextpos: 0].first;
			curdepth--;
		}
	}

	template<typename data_t> long unsigned binary_heap<data_t>::get_pos_depth(const long unsigned &position) const
	{
		long unsigned curdepth = 1;
		long unsigned tempind = 1;
		while(tempind <= position)
		{
			tempind = (tempind+1)*2 -1;
			curdepth++;
		}
		return curdepth;
	}

	template<typename data_t> binary_heap<data_t>::binary_heap()
	: depth(0), n_elements(0), data(NULL)
	{
		data = new pair<data_t, long unsigned>[1];
	}

	template<typename data_t> void binary_heap<data_t>::clear()
	{
		if (data)
			delete[] data;
		
		data = new pair<data_t, long unsigned>[1];
		depth = 0;
		n_elements = 0;
		positions.clear();
	}

	template<typename data_t> binary_heap<data_t>::~binary_heap()
	{
		delete[] data;
	}

	template<typename data_t> bool binary_heap<data_t>::test_consistency() const
	{
		bool all_ok = true;
		if (n_elements == 0)
			return all_ok;
		
		unsigned long tempdepth = get_pos_depth(n_elements-1);
		unsigned long offset = twoe(tempdepth-1) -1;
		unsigned long n_elems = n_elements - offset;
		
		while (tempdepth > 1)
		{
			unsigned long lowo = twoe(tempdepth-2) - 1;
			for (unsigned long i = 0; i < n_elems; i++)
			{
				if(data[offset + i] < data[lowo + i/2])
				{
					cerr << "Binary heap ordering violation at depth " << tempdepth << ", position " << i << "." << endl;
					all_ok = false;
				}
			}
			offset = lowo;
			n_elems = twoe(tempdepth-2);
			tempdepth--;
		}
		return all_ok;
	}

	template<typename d_t> ostream& operator<<(ostream &the_stream, const binary_heap<d_t> &bh)
	{
		the_stream << "Depth: " << bh.depth << endl;
		the_stream << "Number of elements: " << bh.n_elements << endl;
		
		// Build up the list of vertical inter-number distances.
		
		vector<long unsigned> idists;
		
		for (long unsigned i = bh.depth; i > 0; --i)
		{
			long unsigned interdist = twoe(i-1);
			idists.push_back(interdist);
		}
		// Now step through the data and in each line check, which branches are to be written.
		
		long unsigned maxv = twoe(bh.depth > 0 ? bh.depth-1 : 0);
		
		for(long unsigned i = 0; i < maxv; ++i)
		{
			for (long unsigned j = 0; j < bh.depth; ++j)
			{
				if (i%idists[j] == 0)
				{
					long unsigned tmp =bh.depth;
					long unsigned t_i = i;
					while (tmp > j+1)
					{
						t_i = t_i/2;
						tmp--;
					}
					long unsigned offset = twoe(j)-1;
					
					if (offset + t_i < bh.n_elements)
						the_stream << scientific << bh.data[offset+t_i].first << "/" << bh.data[offset+t_i].second << "\t";
				}
				else
				{
					the_stream << "\t" << "\t";
				}
			}
			the_stream << endl;
		}
		
		the_stream << "Indices:" << endl;
		for (map<long unsigned, long unsigned>::const_iterator p = bh.positions.begin(); p != bh.positions.end(); ++p)
		{
			the_stream << p->first << " " << p->second << endl;
		}
		return the_stream;
	}
} // End of namespace simtools.

#endif // BINARY_HEAP_HH
