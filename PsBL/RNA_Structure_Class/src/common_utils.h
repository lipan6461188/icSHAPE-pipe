#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H


#include <iostream>
#include <string>
#include "defines.h"
#include <iterator>
#include <vector>
#include <sstream>
#include <cstring>

using namespace std;

//###################### FILE Utilities ##################################
//Determine if the path <directory>/<file> exists (and is not itself a directory)
bool fileExists(const char* const directory, const char* const file);
//Determine if the file specified by filePath exists and is not itself a directory.
bool fileExists(const char* const filePath);
//Determine if the specified directory exists.
bool dirExists(const char* const fullDirPath);
//! Gets the portion of a path that represents the filename.
//! I.e. removes the directory portion of a path. 
//! If removeExtension is true, the file extension is also removed.
string getFileName(const char * const path, bool removeExtension = false);
//! Gets the directory portion of a path. 
//! It is assumed that the path represents a file.
string getDirName(const char * const path);
//###################### STRING Utilities ##################################
// Utility function that creates a string and fills it by sprintf-formatting the arguments. 
string sfmt(const char* const format, ...);

//! Converts control characters in subject into escape sequences.
//! (useful for debugging strings that might contain control characters.)
//! e.g. "Hello\n" becomes "Hello\\n" 
//! i.e. The linefeed is converted into a literal slash and the letter 'n'
void escapeChars(string &subject);

//! Trim leading whitespace from a string (operates in-place on subject)
void trimLeft(string &subject);

//! Trim trailing whitespace from a string (operates in-place on subject)
void trimRight(string &subject);

//! Trim leading and trailing whitespace from a string (operates in-place on subject)
void trim(string &subject);

//! Converts a string to lower-case (operates on subject in-place).
void toLower(string &subject);

//! Converts a string to upper-case (operates on subject in-place).
void toUpper(string &subject);

//! Find a character in subject (a c-string)
//! Returns the 0-based index where the char is found or string::npos if it is not found.
size_t findchr(const char* const subject, const char find);


//! Returns true if a c-string is NULL or empty ("")
//! Equivalent to (cstr==NULL||strlen(cstr)==0)
inline bool is_blank(const char*const cstr) { return cstr==NULL||*cstr=='\0'; } // if the first character is the null-char, '\0', the string length is 0.

//! Make a copy of a c-string (i.e. a char* or char[] etc) 
//! There's a lot of code duplicated around RNAstructure that does this (and should use
//! this function instead).
//! \returns A pointer to a newly allocated buffer that is a copy of the source buffer.
//!          The returned pointer must be deleted later on.
char* copy_cstr(const char* src);

//! Converts a cstring (char*) into a string, making sure not to
//! dereference a NULL pointer.
//! \return A string copy of cstr or an empty string if cstr is NULL.
inline string as_str(const char*const cstr) { return cstr==NULL?"":cstr; }

//! Attempts to convert a string into an int (e.g. "12" becomes 12)
//! returns true if the conversion succeded or false if it failed (e.g. the value was NOT a number).
bool parseInt(const string &subject, int& resultInt);

//! Attempts to convert a string into an int (e.g. "6.02E23" becomes (double)6.02E23)
//! returns true if the conversion succeded or false if it failed (e.g. the value was NOT a number).
bool parseDbl(const string &subject, double& resultDbl);

//! write vector contents to an output stream. used by join and the 'ofstream<<vector' operator.
template<class T>
void join(ostream &out, const vector<T> &v, const char* const delim=",") {
	if(v.size() > 1) copy(v.begin(),v.end()-1, ostream_iterator<T>(out, delim));
	if(!v.empty()) out << v.back();
}
//! join any vector into a string with elements separated by the specified delimiter.
//! (Useful for cout debugging.)
template<class T>
string join(const vector<T> &v, const char* const delim=",") {
	ostringstream oss; join(oss,v,delim); return oss.str();
}
//! Define the ostream << vector operator. (useful for debugging)
//! e.g. cout << "Vector Contents: " << v << endl;
template <typename T>
std::ostream& operator<<(ostream &out, const vector<T> &v) {
  out << '[';
  join(out,v,", ");
  return out << "]";
}

//! Class used to implement a no-op output stream buffer.
class NullBuffer : public std::streambuf { 	public: int overflow(int c) { return c; } };
//! Class used to implement a no-op output stream.
class NullStream : public std::ostream {
	public: NullStream() : std::ostream(&m_sb) {}
	static NullStream Default;
	private: NullBuffer m_sb;
};

//###################### Pointer and Memory Utilities ##################################

//! The auto_delete class provides a very limited implementation of auto_ptr (precursor to unique_ptr etc)
//! It only does one thing -- deletes a pointer when the auto_delete variable goes out of scope.
//! It is not meant to be copied or passed by value etc.
//! One use is to make it easier to return from a function without having to worry about calling delete for 
//! dynamically allocated objects at the end of the function.
//! Another use is to prevent memory analysis tools (e.g. valgrind) from complaining about 
//! file-level variables (aka global variables) that are dynamically allocated singletons that 
//! live permanently until the end of program execution.
//! With the advent of C++11 there are better ways to handle this and when RNAstructure
//! can routinely be built with c++11 (or later), this class can be removed/replaced.
template<typename T,bool IS_ARRAY=false> 
struct auto_delete {
	auto_delete(T *p): _p(p){}
	~auto_delete() { if (_p!=NULL) { if (IS_ARRAY) delete[] _p; else delete _p;} }
	T& operator*() const { return *_p; }
	operator T*() const { return _p; }
	T* get() const { return _p; }
private:
	T*_p;
	auto_delete(const auto_delete &copy){}  // disable the copy constructor
	auto_delete& operator=(const auto_delete &copy){} // disable the assignment operator
};

#endif // COMMON_UTILS_H
