#ifndef __DEBUG_H
#define __DEBUG_H

#define TRACE

#ifdef TRACE
#define at() printf(" AT: %s:%d: %s()\n", __FILE__, __LINE__, __FUNCTION__); fflush(NULL)
#define intro() printf(" IN: %s:%d: %s() <-----\n", __FILE__, __LINE__, __FUNCTION__); fflush(NULL)
#define outro() printf("OUT: %s:%d: %s() ----->\n", __FILE__, __LINE__, __FUNCTION__); fflush(NULL)
#endif


#endif
