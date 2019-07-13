#ifndef STRING_SPLIT_H
#define STRING_SPLIT_H

#include "pan_type.h"

namespace pan{

//using std::string;

// split a string with a char
StringArray split(const string &String, const char &c);
// a higher efficient method to split s a string
void split(const string &String, const char &c, StringArray &Vec);
// split string with multiple seperators, the splited string saved into Vec
void split(const string &String, const char *seperators, StringArray &Vec);
// split a string with blank spaces
void split(const string &String, StringArray &Vec);


// trimming a string
void trim(string &String, const char &c);
void trim(string &String);
//void trim2(string &String, const char &c);

void trim_left(string &String);
void trim_right(string &String);

void trim_left(string &String, const char &c);
void trim_right(string &String, const char &c);

// find any characters in find_char from String and return the pos
string::size_type findany(const string &String, const char *find_char, string::size_type begin=0);


bool startswith(const char *ref_string, const char *head_string);
inline bool startswith(const string &ref_string, const string &head_string){ return startswith(ref_string.c_str(), head_string.c_str()); }

bool endswith(const char *ref_string, const char *tail_string);
inline bool endswith(const string &ref_string, const string &tail_string){ return endswith(ref_string.c_str(), tail_string.c_str()); }


string reverse_string(const string &input_string);


}
#endif // STRING_SPLIT_H