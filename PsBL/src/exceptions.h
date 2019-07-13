
#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <iostream>
#include <exception>

using namespace std;

namespace pan{

class Invalid_Line: public runtime_error
{
public:

	explicit Invalid_Line(const std::string& what_arg, bool fatal=false): runtime_error(what_arg), fatal(fatal) {};
	explicit Invalid_Line(const char* what_arg, bool fatal=false): runtime_error(what_arg), fatal(fatal) {};

	bool fatal;
};

class Bad_IO: public runtime_error
{
public:

	explicit Bad_IO(const std::string& what_arg, bool fatal=false): runtime_error(what_arg), fatal(fatal) {};
	explicit Bad_IO(const char* what_arg, bool fatal=false): runtime_error(what_arg), fatal(fatal) {};

	bool fatal;
};

class Unexpected_Error: public runtime_error
{
public:
	
	explicit Unexpected_Error(const std::string& what_arg, bool fatal=false): runtime_error(what_arg), fatal(fatal) {};
	explicit Unexpected_Error(const char* what_arg, bool fatal=false): runtime_error(what_arg), fatal(fatal) {};

	bool fatal;
};



}
#endif