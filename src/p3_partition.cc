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
#include "p3_partition.h"
#define MAXB  809600

using namespace std;

//OPEN WRITE FILE
p3_partition::p3_partition (const string & out_prefix, bool write_to_files): partition_out_file(NULL), partition_out_index_file(NULL), partition_out_count_file(NULL), partition_file(NULL) {
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
p3_partition::p3_partition (const string &partition_file_path, const string &range): partition_out_file(NULL), partition_out_index_file(NULL), partition_out_count_file(NULL), partition_file(NULL) {
    // extracting range [start,end]
    size_t pos;
    start = stoi (range, &pos);
    if (pos < range.size()) {
        end = stoi (range.substr(pos+1));
    } else {
        end = start;
    }

    // reading the index file of partitions
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
    if ( end > offsets.size()+1 )
        end = offsets.size()+1;

    partition_file = fopen(partition_file_path.c_str(), "rb");
    fseek(partition_file, offsets[start-1], SEEK_SET);
}

//CALL ON EACH CLUSTER
void p3_partition::add_reads(vector<pair<pair<string, string>, pair<pair<int,int>, int> > > cuts, int p_start, int p_end, string p_ref, int pt_id) {
    size_t pos = ftell(partition_out_file);
    fwrite(&pos, 1, sizeof(size_t), partition_out_index_file);
    fprintf(partition_out_file, "%d %lu %d %d %s %d\n", partition_id, cuts.size(), p_start, p_end, p_ref.c_str(), pt_id);
    for (auto &i: cuts) {
        fprintf(partition_out_file, "%s %s %d %d %d\n", i.first.first.c_str(), i.first.second.c_str(), i.second.first.first, i.second.first.second, i.second.second);
    }
    partition_id++;
}

p3_partition::~p3_partition (){
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

int p3_partition::get_start ()
{
    return p_start;
}

int p3_partition::get_end ()
{
    return p_end;
}

string p3_partition::get_reference ()
{
    return p_ref;
}

int p3_partition::get_id () {
    return partition_id;
}

int p3_partition::get_old_id () {
    return old_id;
}

//read next partition
classified_cuts p3_partition::read_partition() {
    int cut_sz, i;
    char pref[MAXB];
    char name[MAXB], read[MAXB];

    cuts = classified_cuts();

    if (start >= end)
        return cuts;

    partition_count++;
    start++;
    fscanf(partition_file, "%d %d %d %d %s %d\n", &partition_id, &cut_sz, &p_start, &p_end, pref, &old_id);
    p_ref = pref;
    current_cluster.resize(0);
    current_cluster.reserve(cut_sz);

    for (i = 0; i < cut_sz; i++) {
        fgets(pref, MAXB, partition_file);
        int start_pos, end_pos, type;
        sscanf(pref, "%s %s %d %d %d", name, read, &start_pos, &end_pos, &type);
        if (type == 0) {
            cuts.bimodal_cuts.push_back({{string(name), string(read)}, {start_pos, end_pos}});
            cuts.bimodal = true;
        }
        else if (type == 1)
            cuts.left_cuts.push_back({{string(name), string(read)}, {start_pos, end_pos}});
        else if (type == 2)
            cuts.right_cuts.push_back({{string(name), string(read)}, {start_pos, end_pos}});
        else
            cuts.misc_cuts.push_back({{string(name), string(read)}, {start_pos, end_pos}});
        cuts.size += 1;
        current_cluster.push_back({{string(name), string(read)}, {{start_pos, end_pos}, type}});
    }

    return cuts;
}


