
#include "exceptions.hh"

RuntimeError::RuntimeError(string text)
	: message(text)
{}

RuntimeError::RuntimeError(RuntimeError& oldone)
	: message(oldone.message) 
{}

RuntimeError::RuntimeError(const RuntimeError& oldone): message(oldone.message) {}
	
const char* RuntimeError::what() const throw()
{
	return message.c_str();
}
	
RuntimeError::~RuntimeError() throw()
{
}

string tostr(size_t size)
{
	std::stringstream ss;
	ss << size;
	return ss.str();
}

#ifdef USES_PYTHON
void translate_RE(RuntimeError const &re)
{
	//boost::python::object pythonExceptionInstance(re);
	//PyErr_SetObject(ndRuntimeError)
	PyErr_SetString(PyExc_RuntimeError, re.what());
}
#endif // USES_PYTHON
