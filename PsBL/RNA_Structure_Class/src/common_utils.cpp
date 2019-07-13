#include "common_utils.h"
#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

//#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

//Determine if the path <directory>/<file> exists and (is not itself a directory)
bool fileExists(const char* const directory, const char* const file) {
	if (directory==NULL||file==NULL) return false;
	string filePath(directory);
	return fileExists(filePath.append("/").append(file).c_str());
}
bool fileExists(const char* const fullPath) {
	struct stat info;
	return fullPath != NULL && stat(fullPath, &info)==0 && (info.st_mode & S_IFDIR)==0;
}
bool dirExists(const char* const fullDirPath) {
	struct stat info;
	return fullDirPath != NULL && stat(fullDirPath, &info)==0 && (info.st_mode&S_IFDIR)==S_IFDIR;
}

// Gets the portion of a path that represents the filename.
// I.e. removes the directory portion of a path. 
// If removeExtension is true, the file extension is also removed.
string getFileName(const char * const path, bool removeExtension) {
	string filename(path);
	size_t pos = filename.find_last_of("/\\");
	if (pos != string::npos) filename.erase(0, pos + 1);
	if (removeExtension) {
		pos = filename.rfind('.');
		if (pos != string::npos)
			filename.erase(pos);
	}
	return filename;
}

// Gets the directory portion of a path. 
// It is assumed that the path represents a file.
string getDirName(const char * const path) {
	string dir(path);
	size_t pos = dir.find_last_of("/\\");
	if (pos==string::npos) return "."; // e.g. input is "hello.txt" the directory is "." -- the current directory
	dir.resize(pos); // set the length to just before the final slash
	return dir;
}

// Creates string and fills it by sprintf-formatting the arguments.
string sfmt(const char* const format, ...) {
    va_list args;
	int size = strlen(format) + 256;
	char* buf = new char[size];
    va_start(args, format);
    int req_size = vsnprintf(buf, size, format, args);
    va_end(args);
	if (req_size < 0)
		// negative return values indicate an error. 
		req_size = sprintf(buf, "Error formatting arguments: %d", req_size);
	else if (req_size >= size) {  // indicates required buffer size (not including terminating \0)
		// the formatted arguments could not all fit inside the buffer. So resize the buffer to the exact required amount and repeat.
		delete[] buf;
		size = req_size+1; // +1 for \0
		buf = new char[size];
		va_start(args, format);
		vsnprintf(buf, size, format, args);
		va_end(args);
	}
	string ret = buf; // copy the buffer and release it.
	delete[] buf;
	return ret;
}

// Find a char in a c-string
size_t findchr(const char* const subject, const char find) {
	const char* found = strchr(subject, find); // strchr returns a pointer to the char that was found or NULL (if not found)
	return found==NULL ? string::npos : found - subject;
}

// converts control characters in subject into escape sequences.
// (useful for debugging strings that might contain control characters.)
// e.g. "Hello\n" becomes "Hello\\n" 
// i.e. The linefeed is converted into a literal slash and the letter 'n'
void escapeChars(string &subject) {
	string snew;
	snew.reserve((int)(subject.size() * 1.3));
	char numbuf[5];
	for(string::iterator it=subject.begin(); it != subject.end(); it++) {
		char &c = *it;
		if (c < 32 || c > 126)
			switch(c) {
				case '\n': snew+="\\n"; break; 
				case '\r': snew+="\\r"; break;
				case '\t': snew+="\\t"; break;
				case '\0': snew+="\\0"; break;
				default:
					snew += "\\x";
					sprintf(numbuf, "%02X", c);
					snew += numbuf;
					break;
			}
		else
			snew += c;
	}
	subject.swap(snew);
}

void trimLeft(string &s) {
	string::iterator it;
	for (it = s.begin(); it!=s.end(); it++) if (!::isspace(*it)) break; // we are guaranteed to find a non-whitespace character (the final '\0' at s.end() if nothing else)
	s.erase(s.begin(), it);
}
void trimRight(string &s) {
	string::iterator it;
	for (it = s.end()-1; it>=s.begin(); it--) if (!::isspace(*it)) { it++; break; } // start with s.end()-1 to skip the final '\0' at s.end()
	if (it < s.begin()) it++; // increment so we don't delete the NON-whitespace character that was found.
	s.erase(it, s.end());
}
void trim(string &s) {
	trimLeft(s);
	if (!s.empty())
		trimRight(s);
}
// Converts a string to lower-case (operates on subject in-place).
void toLower(string &s) {
	if (s.empty()) return;
	string::iterator it, ite;
	for (it=s.begin(),ite=s.end();it!=ite;it++) 
		*it=tolower(*it);
}
// Converts a string to upper-case (operates on subject in-place).
void toUpper(string &s) {
	if (s.empty()) return;
	string::iterator it, ite;
	for (it=s.begin(),ite=s.end();it!=ite;it++) 
		*it=toupper(*it);
}

 char* copy_cstr(const char* src) {
    if (src == NULL) return NULL;
    char* ptr = new char[strlen(src)+1];
    strcpy(ptr, src);
    return ptr;
}

// Attempts to convert a string into an int (e.g. "12" becomes 12)
// returns true if the conversion succeded or false if it failed (e.g. the value was NOT a number).
bool parseInt(const string& s, int& intVal) {
	stringstream ss(s);
	return !(ss >> intVal).fail();
}
// Attempts to convert a string into an int (e.g. "6.02E23" becomes (double)6.02E23)
// returns true if the conversion succeded or false if it failed (e.g. the value was NOT a number).
bool parseDbl(const string& s, double& dblVal) {
	stringstream ss(s);
	return !(ss >> dblVal).fail();
}

NullStream NullStream::Default; // construct default member.