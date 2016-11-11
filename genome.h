#ifndef __GENOME__
#define __GENOME__

#include<iostream>

using namespace std;

class genome
{
private:
	ifstream fin;
	string reference_name;
	string reference;
//	void load_next(void);

public:
	genome(string);
	~genome();
	string extract(const string&, int &, int &);
	void load_next(void);
	int get_size();
	string get_name();
};

#endif
