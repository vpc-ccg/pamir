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

/**
 * Logger Initiliazer
 * @param name Optional location of Log file
 */
void log_init(string name = "") {
    buffer = new char[20 * 1024 * 1024];
    buffer[0] = '\0';
    buffer_size = 0;
    if (name != "")
        output = fopen(name.c_str(), "w");
    else
        output = stdout;
}

/**
 * Logger closer
 */
void log_close() {
    buffer[buffer_size] = '\0';
    fwrite(buffer, sizeof(char), buffer_size, output);
    if (output != stdout)
        fclose(output);
    delete buffer;
}

/**
 * Logger command to log the given input
 * @param format
 * @param ...
 */
void log(const char *format, ...) {
    va_list args;
    va_start(args, format);
    buffer_size += vsprintf(buffer + buffer_size, format, args);
    va_end(args);
    if (buffer_size > 5 * 1024 * 1024) {
        fwrite(buffer, sizeof(char), buffer_size, output);
        fflush(output);
        buffer_size = 0;
    }
}


