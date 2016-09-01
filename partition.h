#ifndef __PARTITION__
#define __PARTITION__

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <zlib.h>
#include <unordered_set>
using namespace std;

class genome_partition
{
	int distance;
	vector<pair<pair<string, string>, pair<int,int>>> comp;
	unordered_set<string> read_cache;
	unordered_map<string, string> oea_mate;
	int p_start;
	int p_end;
	string p_ref;

	gzFile fp;
	char prev_string[2000];

private:
	bool add_read (string, int, int);
	int fc;

public:
	genome_partition (void);
	genome_partition (const string&, int, const unordered_map<string, string>&);
	~genome_partition (void);
	
	vector<pair<pair<string, string>, pair<int,int>>> get_next (void);
	vector<pair<pair<string, string>, pair<int,int>>> read_partition (const string&, const string&);


	bool has_next (void);
	size_t dump (const vector<pair<pair<string, string>, pair<int, int>>>&, FILE*, int);

	int get_start (void);
	int get_end (void);
	string get_reference (void);
	int output_partition (const string &, const string &);
	int get_cluster_id ();
};


#endif
