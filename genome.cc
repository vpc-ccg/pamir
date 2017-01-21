#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "genome.h"
#include "common.h"
using namespace std;

genome::genome(string filename)
{
	fin.open(filename.c_str());
	reference.reserve(300000000);
	reference_name = "";
	char ch;
	fin.get(ch);
}

genome::~genome()
{
}
void genome::load_all(map<string,string> &chroms)
{
	load_next();
	chroms.emplace(make_pair(reference_name,reference));	
	while(reference!="")
	{
		load_next();
		chroms.emplace(make_pair(reference_name,reference));	
	}
}
void genome::load_next(void)
{	
	reference.clear();
	if (fin.eof())
		return;
	getline(fin, reference_name);
	if(reference_name!="")
		reference_name = reference_name.substr(0, reference_name.find(" "));
	string tmp;
	char ch;
	fin.get(ch);
	while (ch != '>' && !fin.eof())
	{
		reference+=ch;
		getline(fin, tmp);
		reference+=tmp;
		fin.get(ch);
	}

	transform(reference.begin(), reference.end(), reference.begin(), ::toupper);
}

string genome::extract(const string &rname, int &start, int &end)
{
	while (rname != reference_name) {
		load_next();
		
	}
	start = max22(start, 1);
	end = max22(end, 1);
	end = min2(end, (int)reference.size());
	return reference.substr(start-1, end-start+1);
}
char genome::getchar(const string &rname, int &loc)
{
	while (rname != reference_name) {
		load_next();
	}
	return reference[loc-1];
}
int genome::get_size()
{
	return (int)reference.size();
}
string genome::get_name()
{
	return reference_name;
}

