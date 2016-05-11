#ifndef _TIMER_H
#define _TIMER_H

#include <sys/time.h>

class Timer{
public:
	Timer();
	~Timer();

	void start();
	void stop();
	void reset();
	double totaltime();

private:
	double elapsed(const timeval & start, const timeval & stop);

	timeval m_begin, m_end;
	double m_tot;
};

#endif