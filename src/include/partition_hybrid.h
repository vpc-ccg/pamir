#ifndef __PARTITION_HYBRID__
#define __PARTITION_HYBRID__

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <zlib.h>
#include <unordered_set>

#include "sam_processing.h"

using namespace std;

const int MIN_Q = 30;
const int MIN_FREQ = 10;
const int window = 20000;
const int margin = 5;

struct BreakpointCandidate {
    int support = 0;
    vector<pair<string, string> > reads;
    bool in_window = true;
    string chr;
    int window_support = 0;
    int nbr_stretch = 1;
};

class genome_partition_hybrid {
private:
    FILE* partition_out_file;
    FILE* partition_out_index_file;
    FILE* partition_out_count_file;

    string map_file;

    int partition_count = 0;
    int partition_id = 1;

	int distance;

    int p_start;
    int p_end;
    string p_ref;

    //TODO: more todo
    int start;
    int end;
    FILE* partition_file;

    int total;

    vector<pair<pair<string, string>, pair<int,int> > > current_cluster;

    string line;
    string curr_chr;

private:
	void add_read(string, int, int);

public:
	genome_partition_hybrid(const string&, const string&, bool write_to_file = false);
	genome_partition_hybrid(int, const string&,  const string&, const string&, const string&, const string &);
	~genome_partition_hybrid(void);
	genome_partition_hybrid(const string&, bool);

    void cluster_reads(string map_path, int len);

	vector<pair<pair<string, string>, pair<int,int> > > read_partition();
	
	int get_start(void);
	int get_end(void);
	string get_reference(void);
	void output_partitions();
	int get_id();
	int get_total();
    void dump_cluster(vector<BreakpointCandidate>& candidates, int begin, int end_loc,
                                               int breakpoint, string chr);
    int clean_up(vector<BreakpointCandidate>& locs, int offset, bool, int);
    void dump_cluster(vector<pair<string, string> > reads, int, int, string chr);
    pair<string, int> inline fill(vector<BreakpointCandidate>& locs, int genome_offset, ifstream& fin);
};


#endif
