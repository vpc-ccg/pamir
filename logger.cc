#include <iostream>
#include <string>
#include <cstdio>
#include <cctype>
#include <cstdarg>
#include "logger.h"

using namespace std;

FILE *output;
char *buffer;
int buffer_size;

void log_init(string name = "")
{
	buffer = new char[5*1024*1024];
	buffer[0] = '\0';
	buffer_size = 0;
	if (name != "")
		output = fopen(name.c_str(), "w");
	else
		output = stdout;
}


void log_close()
{
	buffer[buffer_size]='\0';
	fwrite (buffer , sizeof(char), buffer_size, output);
	//fprintf(output, "%s", buffer);
	if (output != stdout)
		fclose(output);
	delete buffer;
}

void log(const char* format, ...)
{
	va_list args;
	va_start (args, format);
	buffer_size += vsprintf(buffer+buffer_size, format, args);
	va_end(args);
	if (buffer_size > 5200000) 
	{
		buffer[buffer_size]='\0';
		fwrite (buffer , sizeof(char), buffer_size, output);
		//fprintf(output, "%s", buffer);
		fflush(output);
		buffer_size=0;
		buffer[0] = '\0';
	}
}


