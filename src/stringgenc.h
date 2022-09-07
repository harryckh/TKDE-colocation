#ifndef STRINGGEC_H_
#define STRINGGEC_H_

#include <string>
#include <cstring>
#include <vector>
#include <sstream>
#include <algorithm>
using namespace std;

vector<string> split(string& str, const char* c);
string int2str(int n);
string float2str(double f);
template<class Type> string tostring(Type a);
int str2int(string& s);
double str2float(string& s);
template<class ToType, class FromType> ToType strconvert(FromType t);
string& strtoupper(string& s);
string strtoupper(string s);
string& strtolower(string& s);
string strtolower(string s);

#endif /* STRINGGEC_H_ */
