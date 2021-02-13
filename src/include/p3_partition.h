#ifndef __PARTITION_P3__
#define __PARTITION_P3__

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <zlib.h>
#include <unordered_set>

#include "sketch.h"

using namespace std;

typedef pair<pair<string, string>, pair<int,int> > read_cut_info;

struct p3_read_s {
    string name;
    string sequence;
    range_s range;
    type_en type;
    orientation_en orientation;
};

typedef struct {
    int size = 0;
    bool bimodal = false;
    type_en cluster_type = LONG_INSERTION;
    vector<p3_read_s> bimodal_cuts;
    vector<p3_read_s> left_cuts;
    vector<p3_read_s> right_cuts;
    vector<p3_read_s> misc_cuts;
    vector<p3_read_s> single_peak_cuts;
} classified_cuts;

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
    int total;
    FILE* partition_file;

    int estimated_insertion = -1;

    classified_cuts cuts;
    vector<pair<pair<string, string>, pair<int,int> > > short_reads;

    char prev_string[2000];

    const int INSSIZE=450;

public:
    p3_partition (const string&, bool write_to_file = false);
    p3_partition (const string&,  const string&);
    ~p3_partition (void);

    void add_reads(vector<pair<pair<string, string>, pair<int,int> > > short_reads, vector<p3_read_s> cuts, int p_start, int p_end, string p_ref, int insertion);
    pair<vector<read_cut_info>, classified_cuts> read_partition ();

    int get_start (void);
    int get_end (void);
    string get_reference (void);
    int get_id ();
    int get_total();
    int get_estimated_insertion();
};


#endif
