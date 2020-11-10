#ifndef __PARTITION_P3__
#define __PARTITION_P3__

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <zlib.h>
#include <unordered_set>
using namespace std;

struct classified_cuts {
    int size = 0;
    bool bimodal = false;
    vector<pair<pair<string, string>, pair<int,int> > > bimodal_cuts;
    vector<pair<pair<string, string>, pair<int,int> > > left_cuts;
    vector<pair<pair<string, string>, pair<int,int> > > right_cuts;
    vector<pair<pair<string, string>, pair<int,int> > > misc_cuts;
};

class p3_partition {
private:
    FILE* partition_out_file;
    FILE* partition_out_index_file;
    FILE* partition_out_count_file;

    int partition_count = 0;
    int partition_id = 1;

    int p_start;
    int p_end;
    string p_ref;

    int start;
    int end;
    FILE* partition_file;

    classified_cuts cuts;
    vector<pair<pair<string, string>, pair<int,int> > > short_reads;

    char prev_string[2000];

    const int INSSIZE=450;

public:
    p3_partition (const string&, bool write_to_file = false);
    p3_partition (const string&,  const string&);
    ~p3_partition (void);

    void add_reads(vector<pair<pair<string, string>, pair<int,int> > > short_reads, vector<pair<pair<string, string>, pair<pair<int,int>, int> > > cuts, int p_start, int p_end, string p_ref);
    pair<vector<pair<pair<string, string>, pair<int,int> > >, classified_cuts> read_partition ();

    int get_start (void);
    int get_end (void);
    string get_reference (void);
    int get_id ();
};


#endif
