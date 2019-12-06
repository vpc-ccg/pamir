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

class genome_partition {
private:
	int distance;

    int fc;
    int p_start;
    int p_end;
    string p_ref;

    //TODO: more todo
    int start;
    int end;
    FILE* partition_file;

    //TODO: properly rename it to current_cluster;
	vector<pair<pair<string, string>, pair<int,int> > > current_cluster;

	unordered_set<string> read_cache;    // no IDEA
	unordered_map<string, string> oea_mate;
	map<string, vector<vector<string> > > map_token;
	map<string, string> map_cont;
	// insert orphan into clusters
	map<string, vector<int> > myset;



	gzFile fp;
	char prev_string[2000];

	const int INSSIZE=450;

private:
	bool add_read (string, int, int);

public:
	genome_partition (const string&, const string&);
	genome_partition (const string&, int);
	~genome_partition (void);




	int load_orphan( const string &orphan_contig, const string &oea2orphan);
	int load_oea_mates (const string &mate_file) ;

	vector<pair<pair<string, string>, pair<int,int> > > get_next (void);
	vector<pair<pair<string, string>, pair<int,int> > > read_partition ();

	bool has_next (void);
	size_t dump (const vector<pair<pair<string, string>, pair<int, int>>>&, FILE*, int);

	int get_start (void);
	int get_end (void);
	string get_reference (void);
	void output_partitions();
	int get_cluster_id ();
};


#endif
