#ifndef __PARTITION_P2__
#define __PARTITION_P2__

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <zlib.h>
#include <unordered_set>

#include "sketch.h"

using namespace std;

struct p2_read_s {
    id_t id;
    range_s range;
    type_en type;
    orientation_en orientation;
};

class p2_partition {
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

    int total;
    int estimated_insertion = -1;

    vector<p2_read_s> cut_candidates;
    vector<pair<pair<string, string>, pair<int,int> > > short_reads;


    char prev_string[2000];

    const int INSSIZE=450;

public:
    p2_partition (const string&, bool write_to_file = false);
    p2_partition (const string&,  const string&);
    ~p2_partition (void);

    void add_cuts(vector<pair<pair<string, string>, pair<int,int> > > short_reads,
                  vector<cut> cuts,
                  int p_start, int p_end, string p_ref, int insertion);
    pair<vector<pair<pair<string, string>, pair<int,int> > >, vector<p2_read_s> >
        read_partition();

    int get_start (void);
    int get_end (void);
    string get_reference (void);
    int get_id ();
    int get_old_id();
    int get_total();
    int get_estimated_insertion();
};


#endif
