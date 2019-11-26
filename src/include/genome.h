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

public:
	genome(string);
	~genome();
	string get_bases_at(const string&, int &, int &);
	char get_base_at(const string&, int &);
	void load_next(void);
	int get_size();
	string get_name();
};

#endif
