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
	vector<pair<pair<string, string>, pair<int,int> > > current_cluster;
	unordered_set<string> read_cache;
	unordered_map<string, string> oea_mate;
	// for orphan information
	map<string, vector<vector<string> > > map_token;
	map<string, string> map_cont;
	// insert orphan into clusters
	map<string, vector<int> > myset;

	int p_start;
	int p_end;
	string p_ref;

	gzFile fp;
	char prev_string[2000];
	
	const int INSSIZE=450;

private:
	bool add_read (string, int, int);
	int fc;


public:
	genome_partition (void);
	genome_partition (const string&, int, const unordered_map<string, string>& );
	~genome_partition (void);

	int load_orphan( const string &orphan_contig, const string &oea2orphan);

	vector<pair<pair<string, string>, pair<int,int> > > get_next (void);
	vector<pair<pair<string, string>, pair<int,int> > > read_partition (const string&, const string&);

	bool has_next (void);
	size_t dump (const vector<pair<pair<string, string>, pair<int, int>>>&, FILE*, int);

	int get_start (void);
	int get_end (void);
	string get_reference (void);
	int output_partition (const string &, const string &);
	int get_cluster_id ();
};


#endif
