#include <math.h>
#include <string>
#include <sstream>


//#include <cmath>
#ifndef UTIL_H
#define UTIL_H

#if ! (__GNUG__ == 3) && ! defined (_GLIBCXX_ALGORITHM) && ! defined (_STL_ALGOBASE_H)
template<class T>
inline T max(T a, T b)
{
  return a > b ? a:b;
}

template<class T>
inline T min(T a, T b)
{
  return a < b ? a:b;
}
#endif

#if ! defined (__GNUG__) && ! defined (__HP_aCC)
template<class T>
inline T abs(T a)
{
  return (a > 0) ? a : -a;
}
#endif

template<class T>
inline T sqr(T a)
{
  return a*a;
}

template<class T>
inline T cube(T a)
{
  return a*a*a;
}

namespace cronos {
	template<class T>
	inline T power(T base, int exponent)
	{
		int expo = (exponent > 0) ? exponent : -exponent;
		double ret = 1.;
		for (int mul = 1; mul <= expo; mul++)
			ret *= base;
		return (expo > 0) ? ret : 1./ret;
	}
}//namespace cronos

inline double norm(double a)
{
  return a*a;
}

inline float norm(float a)
{
  return a*a;
}

template<class T>
inline T sign(T a, T b)
{
  return (b >= 0) ? abs(a) : -abs(a);
}

// This could cause compiler errors if used
// in conjunction with (modern?) cmath implementations
// (conflicting redeclaration)
// The error, however, only occurs using c++11
// Since this function is not used (and obsolete?), it can be removed at some point
//
//template<class T>
//inline T log2(T a)
//{
//  const double inv_log2 = 1/log(2.);
//#if (__GNUG__ == 3)
//  return T(inv_log2*log(double(a)));
//#else
//  return T(inv_log2*log(a));
//#endif
//}

double value(const char*, const char* = 0);
bool value_exists(const char*, const char* = 0);
std::string svalue(const char*, const char* = 0);
bool cronos_isnan(double);

double value(std::string, const char* = 0);
bool value_exists(std::string, const char* = 0);
std::string svalue(std::string, const char* = 0);

template <typename T>
std::string ToString(T val)
{
    std::stringstream stream;
    stream << val;
    return stream.str();
}

template <typename T>
std::string ToString(T val, unsigned int padding, std::string pad="0")
{
    std::string str = ToString(val);
	while ( str.size() < padding ) {
		str = pad + str;
    }
    return str;
}


const std::string WHITESPACE = " \n\r\t\f\v";

std::string ltrim(const std::string& s);
std::string rtrim(const std::string& s);
std::string trim(const std::string& s);

#endif


#ifndef PARAMETER_FILE_READER_H
#include "ParameterFileReader.H"
#endif

#ifndef CRONOS_OSTREAM_H
#include "CronosOstream.H"
#endif


