/*
 *  unixtools.h
 *  Missing UNIX functions on WIN32 (MinGW)
 *
 *  iWaQa model framework 2010-2017
 *
 *  SYSTEM/UTILS
 *
 */ 

//assure that we have gettimeofday and usleep also on Windows
#ifdef _WIN32
#ifndef unixtools_h

//rand_r
#include <pthread.h>

// gettimeofday & struct timezone
#include <windows.h>
#include <winsock2.h>
#include <time.h>

//A utility to an uniform mkdir

#include <dir.h>	//TODO: replace for something newer

int mkdir(const char * name, int accessflags)
{
	return mkdir(name);
}

#define EPOCHFILETIME (116444736000000000LLU)

struct timezone
{
	int  tz_minuteswest; /* minutes W of Greenwich */
	int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval * tv, struct timezone * tz)
{
	//function body from Steven Edwards
	//"http://www.winehq.org/pipermail/wine-devel/2003-June/018082.html"
	FILETIME        ft;
	LARGE_INTEGER   li;
	__int64         t;
	static int      tzflag;

	if(tv){
		GetSystemTimeAsFileTime(&ft);
		li.LowPart  = ft.dwLowDateTime;
		li.HighPart = ft.dwHighDateTime;
		t  = li.QuadPart;       /* In 100-nanosecond intervals */
		t -= EPOCHFILETIME;     /* Offset to the Epoch time */
		t /= 10;                /* In microseconds */
		tv->tv_sec  = (long)(t / 1000000UL);
		tv->tv_usec = (long)(t % 1000000UL);
	}
	if(tz){
		if(!tzflag){
			_tzset();
			tzflag++;
		}
		tz->tz_minuteswest = _timezone / 60;
		tz->tz_dsttime = _daylight;
	}
	return 0;
}

int usleep(int useconds)
{
	Sleep(useconds);
	return 0;
}

//define rand_r
int rand_r(unsigned int *seedp)
{
	return rand();
}

#define unixtools_h

#endif	//unixtools_h

#else	//if not _WIN32

	//not on windows
	#include <sys/time.h>
	#include <unistd.h>

#endif	//ifdef _WIN32
