#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include "common.h"
#include "partition_hybrid.h"
#include "sam_processing.h"
#include "sam_parser.h"
#define MAXB  809600

using namespace std;

const int BUCKET_SIZE = 100;
const int MAXN = 1000000;

genome_partition_hybrid::genome_partition_hybrid (const string & out_prefix, bool write_to_files): partition_out_file(NULL), partition_out_index_file(NULL), partition_out_count_file(NULL) {
    if (write_to_files) {
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

genome_partition_hybrid::genome_partition_hybrid (const string &partition_file_path, const string &range, bool write_to_file):genome_partition_hybrid(range, write_to_file) {
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

//void genome_partition_hybrid::cluster_reads(string map_path) {
//    map_file = map_path;
//
//    int search_s = 0;
//    int s = 0;
//    int n = 0;
//
//    vector<int> breakpoints;
//
//    ifstream fin;
//    fin.open(map_path);
//    for (string line; getline(fin, line);) {
//        n++;
//        sam_record record(line);
//
//        if (record.rname != "chr1" || n > 1000000)
//            break;
//
//        if (record.is_fully_mapped) {
//            cerr << record.qname << endl;
//            if (s > 30) {
//                breakpoints.push_back(search_s);
////                exit(0);
//            }
//            s = 0;
//        }
//        else {
//            int loc = record.pos;
//
//            int head_clip = loc + record.head_clip_range;
////            int tail_clip = record.end_pos - record.tail_clip_range;
//            int tail_clip = 0;
//            int insertion_pos = loc + record.insertion_pos;
//
//            int potential_s = max(max(head_clip, tail_clip), insertion_pos);
//            s++;
//            search_s = potential_s;
//        }
//    }
//
//    fin.close();
//    fin.open(map_file);
//
//    int i = 0;
//    partition_id = 1;
//    vector<pair<string, string> > reads;
//
//    int curr_breakpoint = breakpoints[i];
//    int start_range = max(0, curr_breakpoint - 30);
//    int end_range = curr_breakpoint - 10;
//
//    bool f = false;
//    n = 0;
//
//    for (int j = 0; j < breakpoints.size(); j++)
//        cerr << breakpoints[j] << endl;
//
//    for (string line; getline(fin, line);) {
//        if (n > 1000000)
//            break;
//        n++;
//        sam_record record(line);
//
////        cerr << curr_breakpoint << " | " << record.pos << endl;
//
//        if (!record.is_fully_mapped && record.pos >= start_range && record.pos <= end_range) {
//            reads.push_back({record.qname, record.seq});
//            f = true;
//        }
//        else {
//            if (f) {
//                size_t pos = ftell(partition_out_file);
//                fwrite(&pos, 1, sizeof(size_t), partition_out_index_file);
//                fprintf(partition_out_file, "%d %lu %d %d %s\n", partition_id, reads.size(), curr_breakpoint - 500,
//                        curr_breakpoint + 500, "chr1");
//
//                for (auto &i: reads) {
//                    fprintf(partition_out_file, "%s %s %d %d\n", i.first.c_str(), i.second.c_str(), 1, 1);
//                }
//
//                partition_count++;
//                partition_id++;
//
//                i++;
//                if (i >= breakpoints.size())
//                    break;
//                curr_breakpoint = breakpoints[i];
//                start_range = max(0, curr_breakpoint - 130);
//                end_range = curr_breakpoint - 20;
//
//                f = false;
//                reads.clear();
//            }
//        }
//    }
//}

void genome_partition_hybrid::cluster_reads(string map_path) {
    map_file = map_path;

    map<int, float> cnt_new;
    vector<int> locs;

    ifstream fin;
    fin.open(map_path);
    int j = 0;
    for (string line; getline(fin, line);) {
        sam_record record(line);

//        if (record.rname != "chr1")
//            break;

        if ((!record.is_fully_mapped) && ((record.flag & 2048) != 2048) && (record.mapq > 30)) {
            int insertion_pos = record.insertion_pos;

            if (record.head_clip_range != 0)
                locs.push_back(record.pos);

            if (record.tail_clip_range != 0)
                locs.push_back(record.end_pos);

            if (record.insertion_pos != record.pos)
                locs.push_back(record.insertion_pos);
        }
    }

    sort(locs.begin(), locs.end());

//    cerr << "-----------------" << endl;

    vector<int> breakpoints;

    int s = 1;

    int prv_s = locs[0];
    for (int i = 1; i < locs.size(); i++) {
        cerr << locs[i] << endl;
        if (locs[i] == prv_s)
            s++;

        else {
            if (s > 10) {
                breakpoints.push_back(prv_s);
                cerr << prv_s << " | " << s << endl;
            }
            s = 1;
            prv_s = locs[i];
        }
    }

    cerr << "-----------------" << endl;
    cerr << "locs: " << locs.size() << " | breakpoints: " << breakpoints.size() << endl;
    for (int i = 0; i < breakpoints.size(); i++)
        cerr << breakpoints[i] << endl;


    fin.close();
    fin.open(map_path);
    partition_id = 1;

    int i = 0;
    if (!breakpoints.size())
        return;

    int curr = breakpoints[i++];

    vector<pair<string, string> > reads;
    for (string line; getline(fin, line);) {
        sam_record record(line);

//        if ( || ((record.flag & 2048) == 1) )
//            continue;
        if ((record.is_fully_mapped) || (record.flag & 2048) == 2048 || (record.mapq < 30))
            continue;

        if ((record.pos > curr - 300) && (record.pos < curr + 150)) {
            reads.push_back({record.qname, record.seq});
        }
        else if (record.pos > curr + 150) {
            if (reads.size() > 0) {
                size_t pos = ftell(partition_out_file);
                fwrite(&pos, 1, sizeof(size_t), partition_out_index_file);
                fprintf(partition_out_file, "%d %lu %d %d %s\n", partition_id, reads.size(), curr - 300,
                        curr + 300, "chr21");

                for (auto &it: reads) {
                    fprintf(partition_out_file, "%s %s %d %d\n", it.first.c_str(), it.second.c_str(), 1, 1);
                }

                partition_count++;
                partition_id++;

                reads.clear();
            }

            if (i >= breakpoints.size())
                break;

            while (breakpoints[i] - curr < 5) {
                if (i >= breakpoints.size())
                    break;
                curr = breakpoints[i++];
            }
            curr = breakpoints[i++];
        }
    }
    if (reads.size() > 0) {
        size_t pos = ftell(partition_out_file);
        fwrite(&pos, 1, sizeof(size_t), partition_out_index_file);
        fprintf(partition_out_file, "%d %lu %d %d %s\n", partition_id, reads.size(), curr - 300,
                curr + 300, "chr1");

        for (auto &i: reads) {
            fprintf(partition_out_file, "%s %s %d %d\n", i.first.c_str(), i.second.c_str(), 1, 1);
        }

        partition_count++;
        partition_id++;

        reads.clear();
    }
}

genome_partition_hybrid::~genome_partition_hybrid () {
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


int genome_partition_hybrid::get_start ()
{
	return p_start;
}

int genome_partition_hybrid::get_end ()
{
	return p_end;
}

string genome_partition_hybrid::get_reference ()
{
	return p_ref;
}

int genome_partition_hybrid::get_total() {
    return total;
}

int genome_partition_hybrid::get_id () {
	return partition_id;
}

//read next partition
vector <pair<pair < string, string>, pair<int, int>>> genome_partition_hybrid::read_partition() {

	int sz, i;
	char pref[MAXB];
	char name[MAXB], read[MAXB];
//reset:
	if (start > end)
		return vector < pair < pair < string, string >, pair < int, int >> > ();
	partition_count++;
	start++;
	fscanf(partition_file, "%d %d %d %d %s\n", &partition_id, &sz, &p_start, &p_end, pref);
	p_ref = pref;
	current_cluster.resize(0);
	current_cluster.reserve(sz);

	for (i = 0; i < sz; i++) {
		fgets(pref, MAXB, partition_file);
		int loc, support;
		sscanf(pref, "%s %s %d %d", name, read, &support, &loc);
		if (loc != -1)
		    current_cluster.push_back({{string(name), string(read)},
			    					   {support, loc}});
	}
//	if (current_cluster.size() == 0)
//		goto reset;
	return current_cluster;
}