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

#ifndef GROWING_BINARY_SEARCH
#define GROWING_BINARY_SEARCH

#include <iostream>
#include <cmath>
#include "binary_sum_search.hh"

using std::cerr;
using std::endl;
using std::cout;
using std::max;
using std::abs;

namespace simtools
{
	// is_active is a wrapper class that makes it possible for the
	// growing_binary_search tree to treat empty and occupied positions
	// equally without resorting to in-band communication to mark 
	// empty positions.
	template <typename T> class is_active
	{
		public:
			is_active(const T& inval, bool inact)
			: active(inact), _internal(inval)
			{}
			
			is_active()
			: active(false), _internal(0)
			{}
			
			bool operator == (const is_active<T> &right) const
			{
				return ((_internal == right._internal) and (active == right.active));
			}
			
			bool operator < (const is_active<T> &right) const
			{
				return _internal < right._internal;
			}
			
			bool operator > (const is_active<T> &right) const
			{
				return _internal > right._internal;
			}
			
			bool operator < (const double &right) const
			{
				return _internal < right;
			}
			
			bool operator > (const double &right) const
			{
				return _internal > right;
			}
			
			is_active<T>& operator= (const is_active<T> &input)
			{
				_internal = input._internal;
				active = input.active;
				return *this;
			}
			
			is_active<T>& operator= (const T &input)
			{
				_internal = input;
				if (input == 0.)
					active = false;
				else
					active = true;
				
				return *this;
			}
			
			is_active<T>& operator+= (const is_active<T> &input)
			{
				_internal += input._internal;
				active = active | input.active;
				return *this;
			}
			
			is_active<T>& operator-= (const is_active<T> &input)
			{
				_internal -= input._internal;
				active = active | input.active;
				return *this;
			}
			
			is_active<T> operator+ (const is_active<T> &input) const
			{
				is_active<T> temp = *this;
				temp += input;
				return temp;
			}
			
			is_active<T> operator- (const is_active<T> &input) const
			{
				is_active<T> temp = *this;
				temp -= input;
				return temp;
			}
			
			template <typename d_t> friend ostream& operator<<(ostream& the_stream, const is_active<d_t> &ain);
			
			bool active;
			T _internal;
	};

	// growing_binary_search is a binary search tree that allows the 
	// addition and removal of events. In the case of addition, the new 
	// event is appended to the list of already existing possible events. 
	// If this list is at its size limit, the tree is checked for already 
	// removed items at the beginning. If there are any, all entries are 
	// renumbered so that what before was the first nonempty entry of the 
	// list in now the first entry. 
	// In case there are no empty (i.e. previously destroyed) elements at 
	// the beginning, the size of the data field is increased to the next 
	// power of two.
	// Event indices are counted up and are never changed for any event.
	// TODO: Give the tree the optional ability to shrink again.
	
	template<typename data_t> class growing_binary_search : public binary_sum_search< is_active<data_t> >
	{
		public:
			growing_binary_search();
			growing_binary_search(const vector<is_active<data_t> > &input);
			
			unsigned long get_index(const data_t &randnum) const;
			unsigned long insert(const data_t &input);
			void remove(const unsigned long &index);
			
			template<typename d_t> friend ostream& operator<<(ostream &the_stream, const growing_binary_search<d_t> &bis);
			
			void rebuild(const unsigned long &position, const data_t &value, const bool &active = true);	// Change the value in leaf Nr. position to value value and
											// perform all the neccessary modifications to the tree.
			
			long unsigned get_n_entries() const
			{
				return n_entries;
			}
			
			data_t get_propensity() const
			{
				return (binary_sum_search<is_active<data_t> >::data[0].active) ? binary_sum_search<is_active<data_t> >::data[0]._internal : 0.;
			}
			
			template <typename d_t> friend bool check_growing_internals(const growing_binary_search<d_t> &in_hn, const double &in_prop, const unsigned long &in_offs, const unsigned long &in_max, const long unsigned &in_depth);
			template <typename dt> friend bool activity_consistency(const growing_binary_search<dt> &in_hn);
			void propagate_active(const unsigned long &position);
			void full_propagate_active();
		protected:
			unsigned long offset_index;
			unsigned long max_val;	// The maximum number of leaves the tree allows at its current depth.
			long unsigned n_entries;	// The number of leaves that are filled with data.
	};

	template<typename data_t> unsigned long growing_binary_search<data_t>::get_index(const data_t &randnum) const
	{
		if ((randnum < 0.) or (randnum > binary_sum_search<is_active<data_t> >::data[0]._internal))
			throw "Invalid random number passed to get_index.";
		
		unsigned long index = 0;
		long unsigned level = 1;
		data_t temprand = randnum;
		
		while (level < binary_sum_search<is_active<data_t> >::depth)
		{
			long unsigned relind = 2*index;
			long unsigned leftind = relind + twoe(level)-1;
			
			if ((not binary_sum_search<is_active<data_t> >::data[leftind].active) or (temprand > binary_sum_search<is_active<data_t> >::data[leftind]._internal))
			{
				if (not binary_sum_search<is_active<data_t> >::data[leftind+1].active)
				{
					if (binary_sum_search<is_active<data_t> >::data[leftind].active)
					{
						// Again take care of possible floating point fence post problems.
						index = relind;
						temprand = binary_sum_search<is_active<data_t> >::data[leftind]._internal;
					}
					else
					{
						throw "Growing_binary_search binary search ended up between two inactive entries.";
					}
				}
				else
				{
					index = relind+1;
					temprand -= binary_sum_search<is_active<data_t> >::data[leftind]._internal;
					if (temprand < 0.)
						temprand = 0.; // Do this to take care of possible floating point inaccuracy fence post errors.
				}
			}
			else
			{
				index = relind;
			}
			level++;
		}
		
		return index + offset_index;
	}

	template<typename data_t> growing_binary_search<data_t>::growing_binary_search(const vector<is_active<data_t> > &input)
	: binary_sum_search<is_active<data_t> >(input), offset_index(0), max_val(1), n_entries(0.)
	{
		max_val << (binary_sum_search<is_active<data_t> >::depth -1);
	}

	template<typename data_t> growing_binary_search<data_t>::growing_binary_search()
	: binary_sum_search<is_active<data_t> >(vector<is_active<data_t> >()), offset_index(0), max_val(1), n_entries(0.)
	{
	// 	binary_sum_search<is_active<data_t> >::data[0].active = true;
	}

	template<typename data_t> void growing_binary_search<data_t>::propagate_active(const unsigned long &position)
	{
		long unsigned temp_depth = binary_sum_search<is_active<data_t> >::depth;
		unsigned long tempos = position - offset_index;
		long unsigned offset = twoe(temp_depth-1)-1;
		
		while (temp_depth > 1)
		{
			bool a1 = binary_sum_search<is_active<data_t> >::data[offset + tempos].active;
			bool a2 = binary_sum_search<is_active<data_t> >::data[offset + ((tempos%2)? tempos-1:tempos+1)].active;
			
			temp_depth--;
			offset = (offset+1)/2 -1;
			tempos = tempos/2;
			if ((a1 == false) and (a2 == false))
			{
				binary_sum_search<is_active<data_t> >::data[offset + tempos].active = false;
			}
			else
			{
				binary_sum_search<is_active<data_t> >::data[offset + tempos].active = true;
			}
		}
	}

	template<typename data_t> void growing_binary_search<data_t>::full_propagate_active()
	{
		long unsigned temp_depth = binary_sum_search<is_active<data_t> >::depth;
		while (temp_depth > 1)
		{
			long unsigned offset = twoe(temp_depth-1)-1;
			long unsigned n_entries = twoe(temp_depth-2);
			long unsigned ominusone = n_entries-1;
			
			for(long unsigned i = 0; i < n_entries; i++)
			{
				binary_sum_search<is_active<data_t> >::data[ominusone + i].active = (binary_sum_search<is_active<data_t> >::data[offset + 2*i].active) | (binary_sum_search<is_active<data_t> >::data[offset + 2*i+1].active);
			}
			
			temp_depth--;
		}
	}


	template<typename data_t> void growing_binary_search<data_t>::rebuild(const unsigned long &position, const data_t &value, const bool &active)
	{
		unsigned long corr_pos = position - offset_index;
		
		is_active<data_t> temp(value, active);
		binary_sum_search<is_active<data_t> >::rebuild(corr_pos, temp);
	// 	propagate_active(position);
	}

	template<typename data_t> void growing_binary_search<data_t>::remove(const unsigned long &index)
	{
		unsigned long corr_pos = index - offset_index;
		
		if ((corr_pos < 0) or (corr_pos >= max_val))
		{
			cerr << "Out of bound index requested. Must lie between " << offset_index << " and " << offset_index + max_val << endl;
		}
		
		rebuild(index, data_t(0.), false);
		
		// Under no circumstances may n_entries be lowered by removal at the back of the currently filled block. This causes strange
		// index offsets later!
	// 	if (corr_pos == n_entries - 1)
	// 	{
	// 		cout << "Back removal happened." << endl;
	// 		--n_entries;
	// 	}
	}


	template<typename data_t> unsigned long growing_binary_search<data_t>::insert(const data_t &input)
	{
		// First check if we hit the upper bound of our data buffer.
		
		if(n_entries == max_val)
		{
			unsigned long start = binary_sum_search<is_active<data_t> >::depth == 1 ? 0 : ((1 << (binary_sum_search<is_active<data_t> >::depth -1))-1);
			
			is_active<data_t> *spos = &(binary_sum_search<is_active<data_t> >::data[start]);
			is_active<data_t> *real_start = spos;
			
			unsigned long displacement = 0;
			while (not (real_start->active))
			{
				++real_start;
				displacement++;
			}
			
			if (displacement < max(max_val/10, (long unsigned)(1)))
			{
				// We are operating within close margins of the capacity, let's increase it by a factor of 2.
				is_active<data_t>* temp = binary_sum_search<is_active<data_t> >::data;
				binary_sum_search<is_active<data_t> >::data = new is_active<data_t>[max_val*4-1];
				
				is_active<data_t>* curp = real_start;
				for(long unsigned i = 0; i < max_val-displacement; i++, ++curp)
				{
					binary_sum_search<is_active<data_t> >::data[max_val*2-1+i] = *curp;
				}
				
				for(long unsigned i = max_val-displacement; i < 2*max_val; i++)
				{
					binary_sum_search<is_active<data_t> >::data[max_val*2-1+i] = is_active<data_t>(0.,false);
				}
				
				delete[] temp;
				
				max_val *=2;
				
				binary_sum_search<is_active<data_t> >::depth++;
			}
			else
			{
				// There is enough space to reuse our current slice of memory.
				is_active<data_t>* p = spos;
				is_active<data_t>* q = real_start;
				for(long unsigned i = displacement; i < max_val; i++, ++p, ++q)
				{
					*p = *q;
				}
				
				for(long unsigned i = max_val - displacement; i < max_val; i++, ++p)
				{
					*p = is_active<data_t>(0.,false);
				}
			}
			offset_index += displacement;
			n_entries -= displacement;
			
			binary_sum_search<is_active<data_t> >::full_rebuild();
	// 		full_propagate_active();
		}
		rebuild(offset_index + n_entries, input);
		n_entries++;
		
		return offset_index + n_entries - 1;
	}

	template <typename d_t> ostream& operator<<(ostream& the_stream, const is_active<d_t> &ain)
	{
		the_stream << ain._internal << " " << ain.active;
		return the_stream;
	}

	template <typename d_t> ostream& operator<<(ostream& the_stream, const growing_binary_search<d_t> &bis)
	{
		operator<<(the_stream, dynamic_cast<const binary_sum_search<is_active<d_t> >& >(bis));
		
		the_stream << "Number of Entries: " << bis.n_entries << endl;
		the_stream << "Maximum Number: " << bis.max_val << endl;
		the_stream << "Offset: " << bis.offset_index << endl;
		
		return the_stream;
	}

	template <typename d_t> bool check_growing_internals(const growing_binary_search<d_t> &in_hn, const double &in_prop, const unsigned long &in_offs, const unsigned long &in_max, const long unsigned &in_depth)
	{
		if (abs(in_hn.get_propensity() -in_prop) > 1.e-6*in_prop)
			return false;
		
		if (in_hn.offset_index != in_offs)
			return false;
		
		if (in_hn.max_val != in_max)
			return false;
		
		if (in_hn.depth != in_depth)
			return false;
		
		if (not activity_consistency(in_hn))
			return false;
		
		return true;
	}

	template <typename dt> bool activity_consistency(const growing_binary_search<dt> &in_hn)
	{
		// Also do a check of the disabled entries.
		size_t temp_depth = in_hn.depth;
		bool all_ok = true;
		
		while (temp_depth > 1)
		{
			unsigned long n_entries = twoe(temp_depth-1);
			unsigned long offset = n_entries - 1;
			
			unsigned long n_oneup = twoe(temp_depth-2);
			unsigned long ou_offset = n_oneup -1;
			
			for(unsigned long i = 0; i < n_oneup; i++)
			{
				bool a = in_hn.data[offset + 2*i].active;
				bool b = in_hn.data[offset + 2*i+1].active;
				
				bool c = in_hn.data[ou_offset + i].active;
				
				if (not c)
				{
					if(a or b)
					{
						cerr << "Inconsistent activity in growing search tree at depth " << temp_depth-1 << ", position " << i << ". Inactive entry should be active." << endl;
						all_ok = false;
					}
					
					if (in_hn.data[ou_offset + i]._internal > 0)
					{
						cerr << "Nonzero value found for inactive entry at depth " << temp_depth-1 << ", position " << i << ". Is " << scientific << in_hn.data[ou_offset + i]._internal << "." << endl;
						all_ok = false;
					}
				}
				else
				{
					if((not a) and (not b))
					{
						cerr << "Inconsistent activity in growing search tree at depth " <<temp_depth-1 << ", position " << i << ". Active entry should be inactive." << endl;
						all_ok = false;
					}
				}
			}
			temp_depth--;
		}
		return all_ok;
	}
} // End of namespace simtools.

#endif // GROWING_BINARY_SEARCH