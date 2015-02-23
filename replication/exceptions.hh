#ifndef EXCEPTIONS_HH
#define EXCEPTIONS_HH

#include<exception>
#include<sstream>

using std::exception;
using std::string;

string tostr(size_t size);
string tostr(double val);
string tostr(unsigned int val);

struct RuntimeError: std::exception
{
	public:
		RuntimeError(string text);
		
		RuntimeError(RuntimeError& oldone);
		RuntimeError(const RuntimeError& oldone);
		
		const char* what() const throw();
		
		virtual ~RuntimeError() throw();
		
		string get_message();
		
	private:
		string message;
};

#ifdef USES_PYTHON
void translate_RE(RuntimeError const &re);
#endif // USES_PYTHON

#endif // EXCEPTIONS_HH
