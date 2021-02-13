#include <iostream>
#include <string>
#include <cstdio>
#include <map>
#include <vector>
#include <cstdio>
#include <cassert>
#include <zlib.h>
#include <cstring>
#include <cstdlib>
#include "common.h"
#include "p3_partition.h"
#define MAXB  1000000

using namespace std;

typedef int32_t testdef;

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
    if ( end > offsets.size() )
        end = offsets.size();

    total = end - start + 1;

    partition_file = fopen(partition_file_path.c_str(), "rb");
    fseek(partition_file, offsets[start-1], SEEK_SET);
}

//CALL ON EACH CLUSTER
void p3_partition::add_reads(vector<pair<pair<string, string>, pair<int,int> > > short_reads, vector<p3_read_s> cuts, int p_start, int p_end, string p_ref, int insertion) {
    partition_count++;
    size_t pos = ftell(partition_out_file);
    fwrite(&pos, 1, sizeof(size_t), partition_out_index_file);
    fprintf(partition_out_file, "%d %lu %lu %d %d %s %d\n", partition_id, short_reads.size(), cuts.size(), p_start, p_end, p_ref.c_str(), insertion);
    for (auto &i: short_reads) {
        fprintf(partition_out_file, "%s %s %d %d\n", i.first.first.c_str(), i.first.second.c_str(), i.second.first,
                i.second.second);
    }
    for (auto &i: cuts) {
        fprintf(partition_out_file, "%s %s %d %d %u %u\n", i.name.c_str(), i.sequence.c_str(), i.range.start, i.range.end, i.type, i.orientation);
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

int p3_partition::get_total() {
    return total;
}

int p3_partition::get_estimated_insertion() {
    return estimated_insertion;
}

//TODO: check resize(0) performance vs. clear()
//read next partition
pair<vector<read_cut_info>, classified_cuts> p3_partition::read_partition() {
    int cut_sz, sr_sz, i;
    char* pref = (char *)malloc(MAXB);
    char* name = (char *)malloc(MAXB);
    char* read = (char *)malloc(MAXB);

    short_reads.clear();

    classified_cuts cuts;
    cuts.bimodal_cuts.clear();
    cuts.left_cuts.clear();
    cuts.right_cuts.clear();

    if (start <= end) {
        partition_count++;
        start++;
        fscanf(partition_file, "%d %d %d %d %d %s %d\n", &partition_id, &sr_sz, &cut_sz, &p_start, &p_end, pref, &estimated_insertion);
        p_ref = pref;
        short_reads.resize(0);
        short_reads.reserve(sr_sz);

        for (i = 0; i < sr_sz; i++) {
            fgets(pref, MAXB, partition_file);
            int loc, support;
            sscanf(pref, "%s %s %d %d", name, read, &support, &loc);
            short_reads.push_back({{string(name), string(read)}, {support, loc}});
        }

        i = 0;
        while (i < cut_sz) {
            fgets(pref, MAXB, partition_file);
            int start_pos, end_pos, type, orientation;
            sscanf(pref, "%s %s %d %d %d %d", name, read, &start_pos, &end_pos, &type, &orientation);
            range_s range{static_cast<offset_t>(start_pos), static_cast<offset_t>(end_pos)};
            p3_read_s read_tmp{.name = name, .sequence = read, .range = range, .type = static_cast<type_en>(type), .orientation = static_cast<orientation_en>(orientation)};
            if (read_tmp.type == BIMODAL) {
                cuts.bimodal_cuts.push_back(read_tmp);
                cuts.bimodal = true;
                cuts.cluster_type = BIMODAL;
            }
            else if (read_tmp.type == SINGLE_PEAK) {
                cuts.single_peak_cuts.push_back(read_tmp);
                cuts.cluster_type = SINGLE_PEAK;
            }
            else if (read_tmp.type == PARTIAL_LEFT)
                cuts.left_cuts.push_back(read_tmp);
            else if (read_tmp.type == PARTIAL_RIGHT)
                cuts.right_cuts.push_back(read_tmp);
            else
                cuts.misc_cuts.push_back(read_tmp);
            cuts.size += 1;
            i++;
        }
    }

    free(pref);
    free(name);
    free(read);

    return {short_reads, cuts};
}