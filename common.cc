#include <string>
#include <cstdarg>
#include "common.h"

using namespace std;

string reverse_complement ( const string &str )
{
	const char *revComp = "TBGDEFCHIJKLMNOPQRSAUVWXYZ";

	string x;
	for (int i = str.size() - 1; i >= 0; i--)
		x += revComp[str[i] - 'A'];
	return x;
}

string reverse ( const string &str )
{
	string x;
	for (int i = str.size() - 1; i >= 0; i--)
		x += str[i];
	return x;
}

string S (const char* fmt, ...) {
	char *ptr = 0;
    va_list args;
    va_start(args, fmt);
    vasprintf(&ptr, fmt, args);
    va_end(args);
    string s = ptr;
    free(ptr);
    return s;
}
