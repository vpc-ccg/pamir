#ifndef __COMMON__
#define __COMMON__
#include <string>

#define KB  1024LL
#define MB  KB * 1024LL
#define GB  MB * 1024LL

#define VERSION	0x10
#define MAGIC 	(0x07445A00 | VERSION)

#define ERROR(c,...)\
	fprintf(stderr, c"\n", ##__VA_ARGS__)
#define LOG(c,...)\
	fprintf(stderr, c"\n", ##__VA_ARGS__)
#define MAX_CHAR 1000000

using namespace std;

void wo (FILE *, char *, char *, char *);
char checkNs (char *);
string reverse_complement ( const string & );
string reverse ( const string & );
string S (const char* fmt, ...);
int check_AT_GC(const string &, const double &);
/*************************************************/
inline int max2(int a, int b)
{
	return (a>b)?a:b;
}
/*************************************************/
inline int min2(int a, int b)
{
	return (b>a)?a:b;
}
/*************************************************/
inline int max3(int a, int b, int c)
{
	int z = a;
	if (z < b) z=b;
	if (z < c) z=c;
	return z;
}

#endif
