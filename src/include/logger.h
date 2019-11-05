#ifndef __LOGGER__
#define __LOGGER__
#include<string>

using namespace std;

void log_init(string);
void log_close();
void log(const char*, ...);


#endif
