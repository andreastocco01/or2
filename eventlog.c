#include "eventlog.h"
#include <stdio.h>

struct eventlog eventlog = {
    .filename = NULL,
    .file_handle = NULL,
};

int eventlog_initialize(const char* filename)
{
	FILE* file = fopen(filename, "w");
	if (file != NULL) {
		eventlog.filename = filename;
		eventlog.file_handle = file;
		return 0;
	}
	return -1;
}

int eventlog_close()
{
	fflush(eventlog.file_handle);
	fclose(eventlog.file_handle);
	return 0;
}

int eventlog_logdouble(const char* event, int timeinstant, double value)
{
	// I expect the compiler to remove function calls to this
	// when DEBUG is not set. Verify this. TODO
#ifdef DEBUG
	fprintf(eventlog.file_handle, "%s,%d,%lf\n", event, timeinstant, value);
#endif
	return 0;
}
