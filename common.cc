#include <string>
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

