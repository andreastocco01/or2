#ifndef EVENTLOG_H_
#define EVENTLOG_H_

#include <stdio.h>

/**
 * This is a very basic implementation.
 * It is slow by its nature (everytime a new event is added, it is written on the disk).
 * A better implementation would allocate some memory to be used as a buffer.
 * */

struct eventlog {
	const char* filename;
	FILE* file_handle;
};

int eventlog_initialize(const char* filename);
int eventlog_close();

int eventlog_logdouble(const char* event, int timeinstant, double value);

#endif // EVENTLOG_H_
