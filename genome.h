#ifndef __GENOME__
#define __GENOME__

#include<iostream>
#include<map>
using namespace std;

class genome
{
private:
	ifstream fin;
	string reference_name;
	string reference;
	string reserved_refname;
//	void load_next(void);

public:
	genome(string);
	~genome();
	string extract(const string&, int &, int &);
	char getchar(const string&, int &);
	void load_next(void);
	void load_all(map<string,string> &);
	int get_size();
	string get_name();
};

#endif
