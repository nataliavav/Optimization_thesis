//############################################################################
//
//			FILESYSTEM.H
//
///	11/2005   Ioannis Kampolis
//
//	Contains the basic filesystem handling routines (Platform depended)
//	makedir, changedir
//
//############################################################################
//
//
//
//
#ifndef FILESYSTEM_H
#define FILESYSTEM_H


#include <string>
#include <ios>
#include <iostream>
#include <fstream>
#include <vector>
#include "util.h"
#include <string>
#include <streambuf>


std::string GetPassword();
bool changedir(const std::string&);
bool makedir(const std::string&);
bool getworkdir(std::string&);
bool istream_operational(const std::istream&);
void istream_next_char(std::istream&, const char);
void istream_nl(std::istream&);
void istream_get_tokens(std::istream& , std::vector<std::string>&);
template <class T> inline void bin_write(std::ostream& out, const T v)
	{out.write(reinterpret_cast<const char*>(&v), sizeof(v));};
template <class T> inline void bin_read(std::istream& in, T& v)
	{in.read(reinterpret_cast<char*>(&v), sizeof(v));};
void skipLines(std::istream&, const int);

int SystemCall(const std::string);
std::string FileToString(const std::string name);
std::string StreamToString(std::istream&);

#ifndef WIN32
	bool deletedir(const char *path);
#endif

#endif

