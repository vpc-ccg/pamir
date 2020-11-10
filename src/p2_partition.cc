#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cstdio>
#include <cassert>
#include <zlib.h>
#include <cstring>
#include <cstdlib>
#include "common.h"
#include "p2_partition.h"
#define MAXB  809600

using namespace std;

//OPEN WRITE FILE
p2_partition::p2_partition (const string & out_prefix, bool write_to_files): partition_out_file(NULL), partition_out_index_file(NULL), partition_out_count_file(NULL), partition_file(NULL) {
    if(write_to_files){
        partition_out_file = fopen(out_prefix.c_str(), "wb");
        if (partition_out_file == NULL){
            throw("[Genome Partition] Cannot open the file.");
        }
        partition_out_index_file = fopen((out_prefix + ".idx").c_str(), "wb");
        if (partition_out_index_file == NULL){
            throw("[Genome Partition] Cannot open the file.");
        }
        partition_out_count_file = fopen((out_prefix + ".count").c_str(), "w");
        if (partition_out_count_file == NULL){
            throw("[Genome Partition] Cannot open the file.");
        }
    }
}

//READ
p2_partition::p2_partition (const string &partition_file_path, const string &range): partition_out_file(NULL), partition_out_index_file(NULL), partition_out_count_file(NULL), partition_file(NULL) {
    size_t pos;
    start = stoi (range, &pos);
    if (pos < range.size()) {
        end = stoi (range.substr(pos+1));
    } else {
        end = start;
    }

    FILE *partition_file_index = fopen((partition_file_path + ".idx").c_str(), "rb");
    if ( partition_file_index == NULL)
        throw("[Genome Partition] Cannot open the file.");


    Logger::instance().info ("Loading the index file.\n");
    vector<size_t> offsets;
    size_t offset;

    while (fread(&offset, 1, sizeof(size_t), partition_file_index) == sizeof(size_t))
    {
        offsets.push_back(offset);
    }

    fclose(partition_file_index);

    if (start < 1)
        start = 1;
    if ( end > offsets.size() )
        end = offsets.size();

    partition_file = fopen(partition_file_path.c_str(), "rb");
    fseek(partition_file, offsets[start-1], SEEK_SET);
}

//CALL ON EACH CLUSTER
void p2_partition::add_cuts(vector<pair<pair<string, string>, pair<int,int> > > short_reads, vector<pair<string, pair<pair<int, int>, int> > > cuts, int p_start, int p_end, string p_ref, int pt_id) {
    size_t pos = ftell(partition_out_file);
    fwrite(&pos, 1, sizeof(size_t), partition_out_index_file);
    fprintf(partition_out_file, "%d %lu %lu %d %d %s %d\n", partition_id, short_reads.size(), cuts.size(), p_start, p_end, p_ref.c_str(), pt_id);
    for (auto &i: short_reads) {
        fprintf(partition_out_file, "%s %s %d %d\n", i.first.first.c_str(), i.first.second.c_str(), i.second.first,
                i.second.second);
    }
    for (auto &i: cuts) {
        fprintf(partition_out_file, "%s %d %d %d\n", i.first.c_str(), i.second.first.first, i.second.first.second, i.second.second);
    }
    partition_id++;
}

p2_partition::~p2_partition (){
    if (partition_file != NULL){
        fclose(partition_file);
    }
    if (partition_out_file != NULL){
        fclose (partition_out_file);
    }
    if (partition_out_index_file != NULL){
        fclose (partition_out_index_file);
    }
    if (partition_out_count_file != NULL){
        fprintf(partition_out_count_file, "%d\n", partition_count);
        fclose (partition_out_count_file);
    }
}

int p2_partition::get_start ()
{
    return p_start;
}

int p2_partition::get_end ()
{
    return p_end;
}

string p2_partition::get_reference ()
{
    return p_ref;
}

int p2_partition::get_id () {
    return partition_id;
}

int p2_partition::get_old_id () {
    return old_id;
}

//read next partition
vector <pair<string, pair<pair<int, int>, int> > > p2_partition::read_partition() {

    int sr_sz, cut_sz, i;
    char pref[MAXB];
    char name[MAXB], read[MAXB];
    if (start > end)
        return vector <pair<string, pair<pair<int, int>, int> > > ();

    partition_count++;
    start++;
    fscanf(partition_file, "%d %d %d %d %d %s %d\n", &partition_id, &sr_sz, &cut_sz, &p_start, &p_end, pref, &old_id);
    p_ref = pref;
    cut_candidates.resize(0);
    cut_candidates.reserve(cut_sz);

    for (i = 0; i < sr_sz; i++) {
        fgets(pref, MAXB, partition_file);
        int loc, support;
        sscanf(pref, "%s %s %d %d", name, read, &support, &loc);
    }

    for (i = 0; i < cut_sz; i++) {
        fgets(pref, MAXB, partition_file);
        int start_pos, end_pos, type;
        sscanf(pref, "%s %d %d %d", name, &start_pos, &end_pos, &type);
        cut_candidates.push_back({string(name), {{start_pos, end_pos}, type}});
    }

    return cut_candidates;
}