#ifndef TIMEWRAPPER_H
#define TIMEWRAPPER_H

#ifdef _MSC_VER

#define WIN32_LEAN_AND_MEAN

#define NOMINMAX // prevent definition of min() and max() macros, as they clash with min() and max() member functions defined in CronosCode
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64 

// MSVC defines this in winsock2.h!?
struct _timeval {
	long tv_sec;
	long tv_usec;
};

using timeval = struct _timeval;

inline int gettimeofday(timeval* tp, struct timezone* tzp)
{
	// Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
	// This magic number is the number of 100 nanosecond intervals since January 1, 1601 (UTC)
	// until 00:00:00 January 1, 1970 
	static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

	SYSTEMTIME  system_time;
	FILETIME    file_time;
	uint64_t    time;

	GetSystemTime(&system_time);
	SystemTimeToFileTime(&system_time, &file_time);
	time = ((uint64_t)file_time.dwLowDateTime);
	time += ((uint64_t)file_time.dwHighDateTime) << 32;

	tp->tv_sec = (long)((time - EPOCH) / 10000000L);
	tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
	return 0;
}

#else

#include <sys/time.h>

#endif

#endif // TIMEWRAPPER_H