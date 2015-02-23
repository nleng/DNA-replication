// replication - a program to simulate DNA replication in mammalian cells
// 
// Copyright © 2009 - 2012 Daniel Löb <garak@zombiepiratesfromspace.eu>
// 
// replication is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// netdyn is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with replication.  If not, see <http://www.gnu.org/licenses/>.

#ifndef GET_PLOTPOINTS_HH
#define GET_PLOTPOINTS_HH

#include <boost/python.hpp>
#include <vector>
#include <string>

#include <boost/type_traits.hpp>

boost::python::list get_3d_points(boost::python::list in_forks, boost::python::str beadfilename);

#endif // GET_PLOTPOINTS_HH