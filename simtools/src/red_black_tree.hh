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


///////////////////////////////////////////////////////////////////////////////
//
//	This particular red-black tree implementation was heavily inspired by 
//	the english wikipedia page on red-black trees.
//	Wherever possible, the naming conventions from that page were used so
//	that it can serve as a supplementary documentation on how red-black
//	trees work.
//
///////////////////////////////////////////////////////////////////////////////

// Some specialties of this red black tree are:
// - It has iterators
// - Functions find_biggest and find_smallest that return an iterator to the
// 	n-th element from the front or the back. Scalses logarithmically with
// 	number of entries.

#ifndef RED_BLACK_TREE
#define RED_BLACK_TREE

#include <ostream>
#include <list>
#include <limits>
#include <cassert>
#include <bitset>
#include <utility>
#include <string>
#include <stddef.h>

using std::ostream;
using std::endl;
using std::scientific;
using std::list;
using std::bitset;
using std::string;
using std::pair;

namespace simtools
{
	template<typename data_t> struct rb_entry
	{
		rb_entry(const data_t &inval)
		: left(NULL), right(NULL), parent(NULL), me(NULL),
		family_members(1), value(inval), red(true)
		{}
		
		~rb_entry();
		
		bool operator< (const rb_entry<data_t> &right);
		bool operator> (const rb_entry<data_t> &right);
		
		
		bool operator< (const data_t &right);
		bool operator> (const data_t &right);
		
		rb_entry *left;
		rb_entry *right;
		rb_entry *parent;
		rb_entry *me;
		
		unsigned long family_members;	// The number of children + 1
		data_t value;
		bool red;
	};

	template<typename data_t> bool rb_entry<data_t>::operator< (const rb_entry<data_t> &right)
	{
		return value < right.value;
	}

	template<typename data_t> bool rb_entry<data_t>::operator> (const rb_entry<data_t> &right)
	{
		return right.value < value;
	}


	template<typename data_t> bool rb_entry<data_t>::operator< (const data_t &right)
	{
		return value < right;
	}

	template<typename data_t> bool rb_entry<data_t>::operator> (const data_t &right)
	{
		return right < value;
	}


	template<typename data_t> rb_entry<data_t>::~rb_entry()
	{
		if (left)
			delete left;
		if (right)
			delete right;
	}

	template<typename data_t> class red_black_tree
	{
		public:
// 			class const_iterator;
			
			template<typename d_t, bool isconst> class my_iterator_base
			{
				public:
					my_iterator_base(rb_entry<data_t>* input)
					: pointer(input)
					{}
					
#if __GNUC__ >= 4 && __GNUC_MINOR__ >= 5
					explicit operator bool () const
					{
						if (pointer)
							return true;
						return false;
					}
#else
					operator bool () const
					{
						if (pointer)
							return true;
						return false;
					}
#endif
					
					template <bool disconst> bool operator==(const my_iterator_base<d_t, disconst> &right) const
					{
						if (pointer == right.pointer)
							return true;
						return false;
					}
					
					bool operator==(const rb_entry<d_t>* right) const
					{
						if (pointer == right)
							return true;
						return false;
					}
					
					template <bool disconst> bool operator!=(const my_iterator_base<d_t, disconst> &right) const
					{
						if (pointer != right.pointer)
							return true;
						return false;
					}
					
					bool operator!=(const rb_entry<d_t>* right) const
					{
						if (pointer != right)
							return true;
						return false;
					}
					
					friend class red_black_tree<d_t>;
				protected:
					rb_entry<d_t> *pointer;
			};
			
			
			class const_iterator : public my_iterator_base<data_t, true>
			{
				public:
					const_iterator(rb_entry<data_t>* input)
					: my_iterator_base<data_t, true>(input)
					{}
					const data_t* operator->() const
					{
						return &(this->pointer->value);
					}
					
					data_t operator*() const
					{
						return this->pointer->value;
					}
			};
			
			class iterator : public my_iterator_base<data_t, false>
			{
				public:
					iterator(rb_entry<data_t>* input)
					: my_iterator_base<data_t, false>(input)
					{}
					
					data_t* operator->()
					{
						return &(this->pointer->value);
					}
					
					data_t& operator*()
					{
						return this->pointer->value;
					}
					
					operator const_iterator() const
					{
						return const_iterator(this->pointer);
					}
			};
			
			red_black_tree();
			~red_black_tree();
			iterator insert(const data_t &inval);
			template <typename d_t> friend ostream& operator<<(ostream& the_stream, const red_black_tree<d_t> &rbt);
			
			iterator find_biggest(unsigned long N = 0);
			iterator find_smallest(unsigned long N = 0);
			
			const_iterator find_biggest_const(unsigned long N = 0) const;
			const_iterator find_smallest_const(unsigned long N = 0) const;
			const_iterator find_same_const(const data_t &inval) const;
			const_iterator find_smaller_than_const(const data_t &inval) const;
			const_iterator find_bigger_than_const(const data_t &inval) const;
			iterator find_smaller_than(const data_t &inval);
			iterator find_bigger_than(const data_t &inval);
			iterator find_same(const data_t &inval);
			
			void remove(iterator);
			
			unsigned long get_n_elements() const
			{
				return tree ? tree->family_members : 0;
			}
			
			pair<bool, string> check_sortedness() const;
			pair<bool, string> check_family_sizes() const;
			pair<bool, string> check_root_blackness() const;
			pair<bool, string> check_color_correctness() const;
			pair<bool, string> check_black_count() const;
			pair<bool, string> check_backlink_consistency() const;
			
			void clear();
		protected:
			rb_entry<data_t>* grandparent(rb_entry<data_t>* node);
			rb_entry<data_t>* uncle(rb_entry<data_t>* node);
			rb_entry<data_t>* sibling(rb_entry<data_t>* node);
			
			void rotate_left(rb_entry<data_t>* node);
			void rotate_right(rb_entry<data_t>* node);
			void replace_node(rb_entry<data_t>* node, rb_entry<data_t>* replacement);
			
			void insert_case1(rb_entry<data_t>* node);
			void insert_case2(rb_entry<data_t>* node);
			void insert_case3(rb_entry<data_t>* node);
			void insert_case4(rb_entry<data_t>* node);
			void insert_case5(rb_entry<data_t>* node);
			
			void delete_one_child(rb_entry<data_t>* node);
			void delete_case1(rb_entry<data_t>* node);
			void delete_case2(rb_entry<data_t>* node);
			void delete_case3(rb_entry<data_t>* node);
			void delete_case4(rb_entry<data_t>* node);
			void delete_case5(rb_entry<data_t>* node);
			void delete_case6(rb_entry<data_t>* node);
			
			rb_entry<data_t> *tree;
	};

	template<typename data_t> red_black_tree<data_t>::red_black_tree()
		: tree(NULL)
	{
	}

	template<typename data_t> red_black_tree<data_t>::~red_black_tree()
	{
		clear();
	}

	template<typename data_t> void red_black_tree<data_t>::clear()
	{
		delete tree;
		tree = NULL;
	}

	template<typename data_t> typename red_black_tree<data_t>::iterator red_black_tree<data_t>::find_same(const data_t &inval)
	{
		// We're looking for the biggest value entry in the tree that is smaller than inval.
		rb_entry<data_t>* tmp = tree;
		while (tmp != NULL)
		{
			if (tmp->value == inval)
			{
				return iterator(tmp);
			}
			if (tmp->value < inval)
			{
				tmp = tmp->right;
			}
			else
			{
				tmp = tmp->left;
			}
		}
		return iterator(NULL);
	}

	template<typename data_t> typename red_black_tree<data_t>::iterator red_black_tree<data_t>::find_smaller_than(const data_t &inval)
	{
		// We're looking for the biggest value entry in the tree that is smaller than inval.
		rb_entry<data_t>* tmp = tree;
		rb_entry<data_t>* lastfound = NULL;
		while (tmp != NULL)
		{
			if (tmp->value < inval)
			{
				lastfound = tmp;
				tmp = tmp->right;
			}
			else
			{
				tmp = tmp->left;
			}
		}
		if (lastfound)
			if (inval < lastfound->value)
				return iterator(NULL);
		return iterator(lastfound);
	}

	template<typename data_t> typename red_black_tree<data_t>::const_iterator red_black_tree<data_t>::find_smaller_than_const(const data_t &inval) const
	{
		// We're looking for the biggest value entry in the tree that is smaller than inval.
		rb_entry<data_t>* tmp = tree;
		rb_entry<data_t>* lastfound = NULL;
		while (tmp != NULL)
		{
			if (tmp->value < inval)
			{
				lastfound = tmp;
				tmp = tmp->right;
			}
			else
			{
				tmp = tmp->left;
			}
		}
		if (lastfound)
			if (inval < lastfound->value)
				return const_iterator(NULL);
		return const_iterator(lastfound);
	}

	template<typename data_t> typename red_black_tree<data_t>::iterator red_black_tree<data_t>::find_bigger_than(const data_t &inval)
	{
		// We're looking for the smallest value entry in the tree that is bigger than inval.
		rb_entry<data_t>* tmp = tree;
		rb_entry<data_t>* lastfound = NULL;
		while (tmp != NULL)
		{
			if (inval < tmp->value )
			{
				lastfound = tmp;
				tmp = tmp->left;
			}
			else
			{
				tmp = tmp->right;
			}
		}
		if (lastfound)
			if (lastfound->value < inval)
				return iterator(NULL);
		return iterator(lastfound);
	}

	template<typename data_t> typename red_black_tree<data_t>::const_iterator red_black_tree<data_t>::find_bigger_than_const(const data_t &inval) const
	{
		// We're looking for the smallest value entry in the tree that is bigger than inval.
		rb_entry<data_t>* tmp = tree;
		rb_entry<data_t>* lastfound = NULL;
		while (tmp != NULL)
		{
			if (inval < tmp->value )
			{
				lastfound = tmp;
				tmp = tmp->left;
			}
			else
			{
				tmp = tmp->right;
			}
		}
		if (lastfound)
			if (lastfound->value < inval)
				return const_iterator(NULL);
		return const_iterator(lastfound);
	}

	template<typename data_t> typename red_black_tree<data_t>::const_iterator red_black_tree<data_t>::find_same_const(const data_t &inval) const
	{
		// We're looking for the biggest value entry in the tree that is smaller than inval.
		rb_entry<data_t>* tmp = tree;
		while (tmp != NULL)
		{
			if (tmp->value == inval)
			{
				return const_iterator(tmp);
			}
			if (tmp->value < inval)
			{
				tmp = tmp->right;
			}
			else
			{
				tmp = tmp->left;
			}
		}
		return const_iterator(NULL);
	}

	template<typename data_t> typename red_black_tree<data_t>::const_iterator red_black_tree<data_t>::find_biggest_const( unsigned long N) const
	{
		if ((tree == NULL) or (N >= tree->family_members))
		{
			return const_iterator(NULL);
		}
		
		rb_entry<data_t>* tmp = tree;
		unsigned long Nt = N;
		while (tmp != NULL)
		{
			unsigned long right = tmp->right ? tmp->right->family_members: 0;
			if (Nt < right)
			{
				tmp = tmp->right;
				continue;
			}
			else if (Nt == right)
			{
				return const_iterator(tmp);
			}
			else
			{
				Nt -= right + 1;
				tmp = tmp->left;
			}
		}
		return const_iterator(NULL);
	}

	template<typename data_t> typename red_black_tree<data_t>::iterator red_black_tree<data_t>::find_biggest( unsigned long N)
	{
		if ((tree == NULL) or (N >= tree->family_members))
		{
			return iterator(NULL);
		}
		
		rb_entry<data_t>* tmp = tree;
		unsigned long Nt = N;
		while (tmp != NULL)
		{
			unsigned long right = tmp->right ? tmp->right->family_members: 0;
			if (Nt < right)
			{
				tmp = tmp->right;
				continue;
			}
			else if (Nt == right)
			{
				return iterator(tmp);
			}
			else
			{
				Nt -= right + 1;
				tmp = tmp->left;
			}
		}
		return iterator(NULL);
	}

	template<typename data_t> typename red_black_tree<data_t>::iterator red_black_tree<data_t>::find_smallest(unsigned long N)
	{
		if ((tree == NULL) or (N >= tree->family_members))
		{
			return iterator(NULL);
		}
		
		rb_entry<data_t>* tmp = tree;
		unsigned long Nt = N;
		while (tmp != NULL)
		{
			unsigned long left = tmp->left ? tmp->left->family_members: 0;
			if (Nt < left)
			{
				tmp = tmp->left;
				continue;
			}
			else if (Nt == left)
			{
				return iterator(tmp);
			}
			else
			{
				Nt -= left + 1;
				tmp = tmp->right;
			}
		}
		return iterator(NULL);
	}

	template<typename data_t> typename red_black_tree<data_t>::const_iterator red_black_tree<data_t>::find_smallest_const(unsigned long N) const
	{
		if ((tree == NULL) or (N >= tree->family_members))
		{
			return const_iterator(NULL);
		}
		
		rb_entry<data_t>* tmp = tree;
		unsigned long Nt = N;
		while (tmp != NULL)
		{
			unsigned long left = tmp->left ? tmp->left->family_members: 0;
			if (Nt < left)
			{
				tmp = tmp->left;
				continue;
			}
			else if (Nt == left)
			{
				return const_iterator(tmp);
			}
			else
			{
				Nt -= left + 1;
				tmp = tmp->right;
			}
		}
		return const_iterator(NULL);
	}

	template<typename data_t> void red_black_tree<data_t>::remove(iterator tdp)
	{
		// First check, whether or not the node has at most one child. If it doesn't, replace it with one that does.
		
		rb_entry<data_t>* toswap = tdp.pointer;
		rb_entry<data_t>* todelete = tdp.pointer;
		if ((todelete->left) and (todelete->right))
		{
			
			// We have to be really careful here. In order to maintain the validity of all pointers to elements in the
			// tree, we can't operate with simple value swaps. Instead we must make certain, that the memory location
			// that is to be deleted is the one that todelete points to.
			
			rb_entry<data_t>* parent_1 = todelete->parent;
			bool parent_left_1 = false;
			if(parent_1)
			{
				if (parent_1->left == todelete)
					parent_left_1 = true;
			}
		
			rb_entry<data_t>* left_1 = todelete->left;
			rb_entry<data_t>* right_1 = todelete->right;
			unsigned long fmembers_1 = todelete->family_members;
			bool red_1 = todelete->red;
			
			toswap = todelete->right;
			while (toswap->left)
			{
				toswap = toswap->left;
			}
	// 		todelete->value = delnode->value;
			
			bool no_inbetween_element = false;
			if (toswap->parent == todelete)
				no_inbetween_element = true;
			
			if((toswap->parent) and not (no_inbetween_element))
			{
				if (toswap->parent->left == toswap->me)
				{
					toswap->parent->left = todelete;
				}
				else
				{
					toswap->parent->right = todelete;
				}
			}
			todelete->parent = toswap->parent;
			todelete->left = toswap->left;
			if (toswap->left)
			{
				toswap->left->parent = todelete;
			}
			todelete->right = toswap->right;
			if (toswap->right)
			{
				toswap->right->parent = todelete;
			}
			todelete->family_members = toswap->family_members;
			todelete->red = toswap->red;
			
			toswap->parent = parent_1;
			if (parent_1)
			{
				if (parent_left_1)
				{
					parent_1->left = toswap;
				}
				else
				{
					parent_1->right = toswap;
				}
			}
			else
			{
				tree = toswap;
			}
			toswap->left = left_1;
			if (toswap->left)
			{
				toswap->left->parent = toswap;
			}
			toswap->right = right_1;
			if (toswap->right)
			{
				toswap->right->parent = toswap;
			}
			toswap->family_members = fmembers_1;
			toswap->red = red_1;
			
			if (no_inbetween_element)
			{
				if (toswap->right == toswap)
					toswap->right = todelete;
				if (toswap->left == toswap)
					toswap->left = todelete;
				if (toswap->parent == toswap)
					toswap->parent = parent_1;
				if (todelete->parent == todelete)
					todelete->parent = toswap;
			}
		}

		delete_one_child(todelete);
	}

	template<typename data_t> pair<bool, string> red_black_tree<data_t>::check_black_count() const
	{
		if (tree == NULL)
			return pair<bool, string>(true, string("No entries."));
		else
		{
			// Go through the tree and identify all those nodes that have leaves attached.
			
			list<rb_entry<data_t>*> tocheck;
			tocheck.push_back(tree);
			
			list<rb_entry<data_t>*> ntleaves;
				
			typename list<rb_entry<data_t>*>::iterator tmpit = tocheck.begin();
			
			while(true)
			{	
				rb_entry<data_t> *lval = (*tmpit)->left;
				rb_entry<data_t> *rval = (*tmpit)->right;
				
				if (lval)
					tocheck.push_back(lval);
					
				if (rval)
					tocheck.push_back(rval);
				
				if ((not lval) or (not rval))
					ntleaves.push_back(*tmpit);
					
				tocheck.erase(tocheck.begin());
				tmpit = tocheck.begin();
				if (tmpit == tocheck.end())
					break;
			}
			
			bool firstone = true;
			size_t correct_depth = 0;
			for(typename list<rb_entry<data_t>*>::iterator p = ntleaves.begin(); p != ntleaves.end(); ++p)
			{
				rb_entry<data_t>* tmp = (*p);
				size_t counter = 0;
				while (tmp != NULL)
				{
					if (not (tmp->red))
						counter++;
					tmp = tmp->parent;
				}
				if (firstone)
				{
					correct_depth = counter;
					firstone = false;
				}
				else
				{
					if (counter != correct_depth)
					{
						return pair<bool, string>(false, string("Black node counts are inconsistent."));
					}
				}
			}
			
			return pair<bool, string>(true, string("Black node counts are consistent."));
		}
	}

	template<typename data_t> pair<bool, string> red_black_tree<data_t>::check_color_correctness() const
	{
		if (tree == NULL)
			return pair<bool, string>(true, string("No entries."));
		else
		{
			list<rb_entry<data_t>*> tocheck;
			tocheck.push_back(tree);
				
			typename list<rb_entry<data_t>*>::iterator tmpit = tocheck.begin();
			
			while(true)
			{
				rb_entry<data_t> *lval = (*tmpit)->left;
				rb_entry<data_t> *rval = (*tmpit)->right;
				
				if ((*tmpit)->red)
				{
					if ( lval and (lval->red))
					{
						return pair<bool, string>(false, string("Left child of red node is red."));
					}
					
					if (rval and (rval->red))
					{
						return pair<bool, string>(false, string("Right child of right node is red."));
					}
				}
				
				if (lval)
					tocheck.push_back(lval);
					
				if (rval)
					tocheck.push_back(rval);
					
				tocheck.erase(tocheck.begin());
				tmpit = tocheck.begin();
				if (tmpit == tocheck.end())
					break;
			}
			return pair<bool, string>(true, string("Tree colors are consistent."));
		}
	}

	template<typename data_t> pair<bool, string> red_black_tree<data_t>::check_root_blackness() const
	{
		if (tree == NULL)
			return pair<bool, string>(true, string("No entries."));
		else
		{
			if (tree->red)
				return pair<bool, string>(false, string("Root is red."));
		}
		return pair<bool, string>(true, string("Root is black."));
	}

	template<typename data_t> pair<bool, string> red_black_tree<data_t>::check_sortedness() const
	{
		if (tree == NULL)
			return pair<bool, string>(true, string("No entries."));
		else
		{
			list<rb_entry<data_t>*> tocheck;
			tocheck.push_back(tree);
				
			typename list<rb_entry<data_t>*>::iterator tmpit = tocheck.begin();
			
			while(true)
			{
				if (((*tmpit)->family_members == 1) and ( (*tmpit)->left or (*tmpit)->right))
				{
					return pair<bool, string>(true, string("barf."));
				}
				rb_entry<data_t> *lval = (*tmpit)->left;
				rb_entry<data_t> *rval = (*tmpit)->right;
				
				if ( lval and not (lval->value < (*tmpit)->value))
				{
					return pair<bool, string>(false, string("Sorting violation detected in left value."));
				}
				
				if (rval and (rval->value < (*tmpit)->value))
				{
					return pair<bool, string>(false, string("Sorting violation detected in right value."));
				}
				
				if (lval)
					tocheck.push_back(lval);
					
				if (rval)
					tocheck.push_back(rval);
					
				tocheck.erase(tocheck.begin());
				tmpit = tocheck.begin();
				if (tmpit == tocheck.end())
					break;
			}
			return pair<bool, string>(true, string("Tree is sorted."));
		}
	}

	template<typename data_t> pair<bool, string> red_black_tree<data_t>::check_family_sizes() const
	{
		if (tree == NULL)
			return pair<bool, string>(true, string("No entries."));
		else
		{
			list<rb_entry<data_t>*> tocheck;
			tocheck.push_back(tree);
				
			typename list<rb_entry<data_t>*>::iterator tmpit = tocheck.begin();
			
			while(true)
			{
				rb_entry<data_t> *lval = (*tmpit)->left;
				rb_entry<data_t> *rval = (*tmpit)->right;
				
				unsigned long lvm = 0;
				unsigned long rvm = 0;
				
				if (lval)
				{
					lvm = lval->family_members;
				}
				
				if (rval)
				{
					rvm = rval->family_members;
				}
				
				if ((*tmpit)->family_members != rvm + lvm + 1)
				{
					return pair<bool, string>(false, string("Family member inconsistency detected."));
				}
				
				if (lval)
					tocheck.push_back(lval);
					
				if (rval)
					tocheck.push_back(rval);
					
				tocheck.erase(tocheck.begin());
				tmpit = tocheck.begin();
				if (tmpit == tocheck.end())
					break;
			}
			return pair<bool, string>(true, string("Family member count is consistent."));
		}
	}

	template<typename data_t> pair<bool, string> red_black_tree<data_t>::check_backlink_consistency() const
	{
		if (tree == NULL)
			return pair<bool, string>(true, string("No entries."));
		else
		{
			list<rb_entry<data_t>*> tocheck;
			tocheck.push_back(tree);
			
			if (tree->parent != NULL)
				return pair<bool, string>(false, string("Tree root parent exists. This is wrong."));
				
			typename list<rb_entry<data_t>*>::iterator tmpit = tocheck.begin();
			
			while(true)
			{
				rb_entry<data_t> *lval = (*tmpit)->left;
				rb_entry<data_t> *rval = (*tmpit)->right;
				
				if (lval)
				{
					if (lval->parent != (*tmpit))
						return pair<bool, string>(false, string("Left link parent does not point to actual parent."));
				}
				
				if (rval)
				{
					if (rval->parent != (*tmpit))
						return pair<bool, string>(false, string("Right link parent does not point to actual parent."));
					
				}
				
				if (lval)
					tocheck.push_back(lval);
					
				if (rval)
					tocheck.push_back(rval);
					
				tocheck.erase(tocheck.begin());
				tmpit = tocheck.begin();
				if (tmpit == tocheck.end())
					break;
			}
		}
		return pair<bool, string>(true, string("Parent links are all consistent."));
	}

	template<typename data_t> void red_black_tree<data_t>::rotate_left(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* newroot = node->right;
		newroot->parent = node->parent;
		
		node->right = newroot->left;
		if (node->right)
			node->right->parent = node;
		
		newroot->left = node;
		newroot->left->parent = newroot;
		
		if (newroot->parent)
		{
			if (node == newroot->parent->left)
			{
				newroot->parent->left = newroot;
			}
			else
			{
				newroot->parent->right = newroot;
			}
		}
		else
			tree = newroot;
		
		node->family_members -= newroot->family_members;
		if (node->right)
		{
			newroot->family_members -= node->right->family_members;
			node->family_members += node->right->family_members;
		}
		newroot->family_members += node->family_members;
	}

	template<typename data_t> void red_black_tree<data_t>::rotate_right(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* newroot = node->left;
		newroot->parent = node->parent;
		
		node->left = newroot->right;
		if (node->left)
			node->left->parent = node;
		
		newroot->right = node;
		newroot->right->parent = newroot;
		
		if (newroot->parent)
		{
			if (node == newroot->parent->left)
			{
				newroot->parent->left = newroot;
			}
			else
			{
				newroot->parent->right = newroot;
			}
		}
		else
			tree = newroot;
		
		node->family_members -= newroot->family_members;
		if (node->left)
		{
			newroot->family_members -= node->left->family_members;
			node->family_members += node->left->family_members;
		}
		newroot->family_members += node->family_members;
	}

	template<typename data_t> void red_black_tree<data_t>::insert_case1(rb_entry<data_t>* node)
	{
		if(node->parent == NULL)
			node->red = false;
		else
			insert_case2(node);
	}

	template<typename data_t> void red_black_tree<data_t>::insert_case2(rb_entry<data_t>* node)
	{
		if(not(node->parent->red))
			return;
		else
			insert_case3(node);
	}

	template<typename data_t> void red_black_tree<data_t>::insert_case3(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* u = uncle(node);
		rb_entry<data_t>* gp;
		if((u != NULL) and (u->red))
		{
			node->parent->red = false;
			u->red = false;
			gp = grandparent(node);
			gp->red = true;
			insert_case1(gp);
		}
		else
		{
			insert_case4(node);
		}
	}

	template<typename data_t> void red_black_tree<data_t>::insert_case4(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* gp = grandparent(node);
		
		if ((node == node->parent->right) and (node->parent == gp->left))
		{
			rotate_left(node->parent);
			node = node->left;
		}
		else if ((node == node->parent->left) and (node->parent == gp->right))
		{
			rotate_right(node->parent);
			node = node->right;
		}
		insert_case5(node);
	}

	template<typename data_t> void red_black_tree<data_t>::insert_case5(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* gp = grandparent(node);
		
		node->parent->red = false;
		gp->red = true;
		if ((node == node->parent->left) and (node->parent == gp->left))
		{
			rotate_right(gp);
		}
		else
		{
			rotate_left(gp);
		}
	}

	template<typename data_t> void red_black_tree<data_t>::delete_one_child(rb_entry<data_t>* node)
	{
		// At this point, we know that node has at most one child.
		rb_entry<data_t>* child = (node->right == NULL) ? node->left : node->right;
		
		if (child == NULL)
		{
			if (not (node->red))
				delete_case1(node);
			replace_node(node, child);
		}
		else
		{
			replace_node(node, child);
			if(not (node->red))
			{
				if (child->red)
				{
					child->red = false;
				}
				else
				{
					delete_case1(child);
				}
			}
		}
		node->left = NULL;
		node->right = NULL;
		delete node;
	}

	template<typename data_t> void red_black_tree<data_t>::delete_case1(rb_entry<data_t>* node)
	{
		if (node->parent != NULL)
			delete_case2(node);
	}

	template<typename data_t> void red_black_tree<data_t>::delete_case2(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* sib = sibling(node);
		
		if (sib->red)
		{
			node->parent->red = true;
			sib->red = false;
			if (node == node->parent->left)
				rotate_left(node->parent);
			else
				rotate_right(node->parent);
		}
		delete_case3(node);
	}

	template<typename data_t> void red_black_tree<data_t>::delete_case3(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* sib = sibling(node);
		
		if ((not node->parent->red) and (not sib->red) and ((not sib->left) or (not sib->left->red)) and ((not sib->right) or (not sib->right->red)))
		{
			sib->red = true;
			delete_case1(node->parent);
		}
		else
		{
			delete_case4(node);
		}
	}

	template<typename data_t> void red_black_tree<data_t>::delete_case4(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* sib = sibling(node);
		
		if ((node->parent->red) and (not sib->red) and ((not sib->left) or (not sib->left->red)) and ((not sib->right) or (not sib->right->red)))
		{
			sib->red = true;
			node->parent->red = false;
		}
		else
		{
			delete_case5(node);
		}
	}

	template<typename data_t> void red_black_tree<data_t>::delete_case5(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* sib = sibling(node);
		
		if (not sib->red)
		{
			if ((node == node->parent->left) and ((not sib->right) or (not sib->right->red)) and ( (sib->left) and (sib->left->red)))
			{
				sib->red = true;
				sib->left->red = false;
				rotate_right(sib);
			}
			else if ((node == node->parent->right) and ((not sib->left) or (not sib->left->red)) and ((sib->right) and (sib->right->red)))
			{
				sib->red = true;
				sib->right->red = false;
				rotate_left(sib);
			}
		}
		delete_case6(node);
	}

	template<typename data_t> void red_black_tree<data_t>::delete_case6(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* sib = sibling(node);
		
		sib->red = node->parent->red;
		node->parent->red = false;
		
		if (node == node->parent->left)
		{
			sib->right->red = false;
			rotate_left(node->parent);
		}
		else
		{
			sib->left->red = false;
			rotate_right(node->parent);
		}
	}

	template<typename data_t> void red_black_tree<data_t>::replace_node(rb_entry<data_t>* node, rb_entry<data_t>* replacement)
	{
		if (replacement != NULL)
		{
			replacement->parent = node->parent;
		}
		
		rb_entry<data_t>* tmp = node->parent;
		while(tmp != NULL)
		{
			tmp->family_members--;
			tmp = tmp->parent;
		}
		if (node->parent)
		{
			if (node == node->parent->left)
			{
				node->parent->left = replacement;
			}
			else
			{
				node->parent->right = replacement;
			}
		}
		else
		{
			tree = replacement;
		}
	}

	template<typename data_t> rb_entry<data_t>* red_black_tree<data_t>::grandparent(rb_entry<data_t>* node)
	{
		if((node != NULL) and (node->parent != NULL))
			return node->parent->parent;
		else
			return NULL;
	}

	template<typename data_t> rb_entry<data_t>* red_black_tree<data_t>::uncle(rb_entry<data_t>* node)
	{
		rb_entry<data_t>* gp = grandparent(node);
				
		if(gp == NULL)
			return NULL;
		if (node->parent == gp->left)
			return gp->right;
		else
			return gp->left;
	}

	template<typename data_t> rb_entry<data_t>* red_black_tree<data_t>::sibling(rb_entry<data_t>* node)
	{
		if (node->parent)
		{
			if (node == node->parent->left)
			{
				return node->parent->right;
			}
			else
			{
				return node->parent->left;
			}
		}
		else
			return NULL;
	}

	template<typename data_t> typename red_black_tree<data_t>::iterator red_black_tree<data_t>::insert(const data_t &inval)
	{
		rb_entry<data_t>* newval = new rb_entry<data_t>(inval);
		newval->me = newval;
		
		rb_entry<data_t>** tptr = &tree;
		rb_entry<data_t>** parent = &tree;
		
		while(true)
		{
			if((*tptr) == NULL)
			{
				newval->parent = *parent;
				*tptr = newval;
				break;
			}
			else
			{
				parent = tptr;
				(*parent)->family_members++;	
				if (newval->value < (*parent)->value)
				{
					tptr = &((*tptr)->left);
				}
				else
				{
					tptr = &((*tptr)->right);
				}
			}
		}
		insert_case1(*tptr);

		return iterator(newval->me);
	}

	template <typename d_t> ostream& operator<<(ostream& the_stream, const red_black_tree<d_t> &rbt)
	{
		if (rbt.tree == NULL)
			the_stream << "No entries." << endl;
		else
		{
			list<rb_entry<d_t>*> toplot;
			toplot.push_back(rbt.tree);
			
			typename list<rb_entry<d_t>*>::iterator tmpit = toplot.begin();

			unsigned long stage_max = 1;
			unsigned long counter = 0;
			while(true)
			{
				if (counter == stage_max)
				{
					the_stream << endl;
					stage_max *= 2;
				}
				the_stream << (*tmpit)->family_members << "\t" << (*tmpit)->red << "\t" << (*tmpit)->value << endl;
				
				rb_entry<d_t> *lval = (*tmpit)->left;
				rb_entry<d_t> *rval = (*tmpit)->right;
				
				if (lval)
					toplot.push_back(lval);
				
				if (rval)
					toplot.push_back(rval);
				
				toplot.erase(toplot.begin());
				tmpit = toplot.begin();
				if (tmpit == toplot.end())
					break;
				
				counter++;
			}
		}
		return the_stream;
	}
} // End of namespace simtools.

#endif // RED_BLACK_TREE
