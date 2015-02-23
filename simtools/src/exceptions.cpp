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

#include "exceptions.hh"

simtools::RuntimeError::RuntimeError(string text)
	: message(text)
{}

simtools::RuntimeError::RuntimeError(RuntimeError& oldone)
	: message(oldone.message) 
{}

simtools::RuntimeError::RuntimeError(const RuntimeError& oldone): message(oldone.message) {}
	
const char* simtools::RuntimeError::what() const throw()
{
	return message.c_str();
}
	
simtools::RuntimeError::~RuntimeError() throw()
{
}

#ifdef USES_PYTHON
// This is necessary to pass through thrown exception messages to python.
void simtools::translate_RE(RuntimeError const &re)
{
	PyErr_SetString(PyExc_RuntimeError, re.what());
}
#endif // USES_PYTHON
