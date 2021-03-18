#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include "common.h"
#include "progressbar.h"
#include "partition_hybrid.h"

#define MAXB  809600

using namespace std;

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

void genome_partition_hybrid::dump_cluster(vector<BreakpointCandidate>& candidates, int begin, int end_loc,
                                           int breakpoint, string chr) {
//    //cerr << "Dumping\t" << begin << "-" << end_loc << endl;
    string cluster_content = "";
    int reads_size = 0;
    for (int i = begin; i <= end_loc; i++) {
//        //cerr << i << " | " << candidates[i].support << " | " << candidates[i].reads.size() << endl;
        for (int j = 0; j < candidates[i].reads.size(); j++) {
            reads_size++;
            cluster_content += candidates[i].reads[j].first + " " + candidates[i].reads[j].second + " 1 1\n";
        }
    }
//    //cerr << "Reads: " << reads_size << endl;
    if (reads_size == 0)
        return;
    size_t pos = ftell(partition_out_file);
    fwrite(&pos, 1, sizeof(size_t), partition_out_index_file);
    fprintf(partition_out_file, "%d %d %d %d %s\n", partition_id, reads_size, breakpoint - 150,
            breakpoint + 150, chr.c_str());
    fprintf(partition_out_file, "%s", cluster_content.c_str());
    partition_count++;
    partition_id++;
}

void genome_partition_hybrid::dump_cluster(vector<pair<string, string> > reads, int range_start, int range_end, string chr) {
    // cerr << "Dumping" << endl;
	if (reads.size() == 0)
		return;
    size_t pos = ftell(partition_out_file);
    fwrite(&pos, 1, sizeof(size_t), partition_out_index_file);
    fprintf(partition_out_file, "%d %d %d %d %s\n", partition_id, reads.size(), range_start,
            range_end, curr_chr.c_str());
    for (int i = 0; i < reads.size(); i++) {
        fprintf(partition_out_file, "%s %s 1 1\n", reads[i].first.c_str(), reads[i].second.c_str());
    }
    partition_count++;
    partition_id++;
}

//void genome_partition_hybrid::clean_up(vector<BreakpointCandidate>& locs, int offset) {
//    int curr_candidate = -1;
//    int curr_support = -1;
//    int prv_breakpoint = -1;
//    vector<pair<string, string> > curr_reads;
//    int sum = 0;
//    for (int i = -5; i < 6; i++) {
//        if (i + margin > 0) {
//            sum += locs[i + margin].support;
//        }
//    }
//    locs[margin].window_support = sum;
//    for (int i = margin + 1; i < window - 5; i++) {
//        locs[i].window_support = locs[i - 1].window_support - locs[i - 6].support + locs[i + 5].support;
//        if (locs[i].window_support > 5)
//            locs[i].nbr_stretch = 0;
//        else
//            locs[i].nbr_stretch = locs[i-1].nbr_stretch + 1;
//    }
////    for (int i = margin; i < window; i++) {
////        //cerr << offset + i << ": " << locs[i].support << " | " << locs[i].window_support << " | " << locs[i].br << endl;
////        if (curr_candidate == -1 || locs[i].chr == locs[curr_candidate].chr) {
////            if (locs[i].support > MIN_FREQ && locs[i].support > curr_support) {
////                curr_candidate = i;
////                curr_support = locs[i].support;
////            }
////        }
////        else {
////            if (curr_candidate != -1 && locs[i].in_window && curr_candidate - prv_breakpoint > 10) {
////                dump_cluster(locs, max(0, curr_candidate - 150), curr_candidate,
////                             offset + curr_candidate, locs[curr_candidate].chr);
////                prv_breakpoint = curr_candidate;
////            }
////            curr_candidate = -1;
////            curr_support = -1;
////        }
////    }
//
//    int last_point = 0;
//    int cnt = 0;
//    for (int i = locs.size() - 5; i > locs.size() - margin; i--) {
//       if (locs[i].nbr_stretch > 10)
//           last_point = i - 10;
//    }
//
//    int stretch_start = -1, stretch_end = -1;
//    vector<pair<string, string> > reads;
//    for (int i = 0; i < last_point; i++) {
////        //cerr << offset + i << ": " << locs[i].support << " | " << locs[i].window_support << " | " << locs[i].br << endl;
//        if (locs[i].nbr_stretch == 0) {
//            if (stretch_start == -1) {
//                stretch_start = i;
//                reads.insert(reads.end(), locs[i].reads.begin(), locs[i].reads.end());
//            }
//            else
//                reads.insert(reads.end(), locs[i].reads.begin(), locs[i].reads.end());
//        }
//        else {
//            if (stretch_start != -1) {
//                stretch_end = i - 1;
//                int breakpoint_location = offset + (stretch_end - stretch_start)/2;
//                dump_cluster(reads, breakpoint_location, locs[stretch_start].chr);
//                stretch_start = -1;
//                stretch_end = -1;
//                reads.clear();
//            }
//        }
//    }
//
//    int j = 0;
//    for (int i = last_point - 10; i < locs.size(); i++) {
//        locs[j] = locs[i];
//        j++;
//        locs[i].in_window = false;
//    }
//    for (int i = j; i < locs.size(); i++) {
//        BreakpointCandidate tmp;
//        locs[i] = tmp;
//    }
//}

int genome_partition_hybrid::clean_up(vector<BreakpointCandidate>& locs, int offset, bool process_full, int next_offset) {
    int curr_candidate = -1;
    int curr_support = -1;
    int prv_breakpoint = -1;
    vector<pair<string, string> > curr_reads;

    int sum = 0;
    for (int i = 0; i < 11; i++) {
        sum += locs[i].support;
    }
    locs[5].window_support = sum;
   if (locs[5].window_support >= 10)
       locs[5].nbr_stretch = 0;

    for (int i = 6; i < window; i++) {
        locs[i].window_support = locs[i - 1].window_support - locs[i - 6].support + locs[i + 5].support;
        if (locs[i].window_support >= 10)
            locs[i].nbr_stretch = 0;
        else
            locs[i].nbr_stretch = locs[i-1].nbr_stretch + 1;
    }

    int last_point = window;
    if (!process_full) {
        int start = window;
		// cerr << "start: " << start << endl; 
        if (next_offset - offset < window) {
            start = next_offset - offset;
		}
		// cerr << "start: " << start << endl;
        int cnt = 0;
		last_point = start - 1;
        for (int i = start - 5; i > 0; i--) {
            if (locs[i].nbr_stretch > 10) {
                last_point = i - 10;
				// cerr << "last_point: " << last_point << endl;
                break;
            }
        }
    }

	// cerr << "end: " << last_point << endl;

    int stretch_start = -1, stretch_end = -1;
    vector<pair<string, string> > reads;
    for (int i = 0; i < last_point; i++) {
    //    cerr << offset + i << ": " << locs[i].support << " | " << locs[i].window_support << " | " << locs[i].nbr_stretch << endl;
        if (locs[i].nbr_stretch == 0) {
            if (stretch_start == -1) {
                stretch_start = i;
                // cerr << "Found: " << stretch_start << endl;
                reads.insert(reads.end(), locs[i].reads.begin(), locs[i].reads.end());
            }
            else {
                reads.insert(reads.end(), locs[i].reads.begin(), locs[i].reads.end());
            }
        }
        else {
            if (stretch_start != -1) {
                stretch_end = i - 1;
                // cerr << stretch_start << " | " << stretch_end << endl;
                int breakpoint_location = offset + (stretch_end + stretch_start)/2;
                dump_cluster(reads, breakpoint_location - 150, breakpoint_location + 150, locs[stretch_start].chr);
				// if (reads.size() > 0)
				// 	cerr << "size: " << stretch_end - stretch_start << endl;
                stretch_start = -1;
                stretch_end = -1;
                reads.clear();
            }
        }
    }

    int j = 0;
    if (!process_full) {
        for (int i = last_point; i < window; i++) {
            locs[j] = locs[i];
			if (j < 10)
            	locs[j].nbr_stretch = 1;
            j++;
            locs[i].in_window = false;
        }
    }

//    //cerr << "j: " << j << endl;
    for (int i = j; i < window; i++) {
        BreakpointCandidate tmp;
        locs[i] = tmp;
    }

    return j;
}

pair<string, int> inline genome_partition_hybrid::fill(vector<BreakpointCandidate>& locs, int genome_offset,
                                         ifstream& fin) {
    do {
        sam_record record(line);

	// && (record.mapq > MIN_Q)
	// && (record.flag & 8 == 8)
        if ((!record.is_fully_mapped) && ((record.flag & 2048) != 2048) && (record.mapq > MIN_Q)) {
			// cerr << record.qname << "\t s: " << record.pos << " | e: " << record.end_pos << " | hc: " << record.head_clip_range << " | tc: " <<
				// record.tail_clip_range <<  " | i: " << record.insertion_pos << "(" << genome_offset << ")" << endl;

			// if (record.head_clip_range != 0 && record.tail_clip_range != 0)
			// 	continue;

            if (record.rname != curr_chr) {
                return {record.rname, record.pos};
            }

            if (record.pos >= genome_offset + window || record.end_pos >= genome_offset + window) {
                return {record.rname, record.pos};
            }

            if (record.head_clip_range != 0) {
				//cerr << "*1" << endl; 
                // if (record.pos - genome_offset < 0)
                //     //cerr << "1* " << record.qname << " | " << genome_offset << record.pos << endl;
                locs[record.pos - genome_offset].support++;
                locs[record.pos - genome_offset].reads.push_back({record.qname, record.seq});
                locs[record.pos - genome_offset].chr = record.rname;
            }

            if (record.tail_clip_range != 0) {
				//cerr << "-1" << endl; 
                // if (record.end_pos - genome_offset - 1 < 0)
                //     //cerr << "2* " << record.qname << " | " << genome_offset << record.pos << endl;
//                if (record.end_pos - genome_offset - 1 >= window)
//                    return {record.rname, record.pos};
                locs[record.end_pos - genome_offset - 1].support++;
				// if (record.end_pos - genome_offset - 1 < 0 || record.end_pos - genome_offset - 1 > window)
				// cerr << "*** " << record.end_pos - genome_offset - 1 << endl;
                locs[record.end_pos - genome_offset - 1].reads.push_back({record.qname, record.seq});
                locs[record.end_pos - genome_offset - 1].chr = record.rname;
            }

            if (record.insertion_pos != record.pos) {
				//cerr << "+1" << endl; 
                // if (record.insertion_pos - genome_offset < 0)
                //     //cerr << "3* " << record.qname << " | " << genome_offset << record.pos << endl;
//                if (record.insertion_pos - genome_offset >= window)
//                    return {record.rname, record.pos};
                locs[record.insertion_pos - genome_offset].support++;
                locs[record.insertion_pos - genome_offset].reads.push_back({record.qname, record.seq});
                locs[record.insertion_pos - genome_offset].chr = record.rname;
            }

//            if (record.pos >= 16858299 && record.pos <= 16858301 || record.end_pos >= 16858299 && record.end_pos <= 16858301 ||
//                    record.insertion_pos >= 16858299 && record.insertion_pos <= 16858301)
//                //cerr << "# " << record.qname << endl;
        }
    }
    while (getline(fin, line));

    return {"end", INT_MAX};
}

void genome_partition_hybrid::cluster_reads(string map_path, int len) {
    map_file = map_path;
	map<string, int> chr_len;

	ProgressBar progress(80);

    ifstream fin;
    fin.open(map_path);

    int genome_offset, array_offset = 5;

    vector<BreakpointCandidate> locs;
    for (int i = 0; i < window + margin; i++) {
        BreakpointCandidate tmp;
        locs.push_back(tmp);
    }

	bool show_progress = true;
    string next_chr;
    int next_genome_offset;
    bool process_full = false;
    string line;
	int length;
	uint64_t sum = 0;
    getline(fin, line);
	while (line[0] == '@') {
		if (line.substr(0, 3) == "@SQ") {
			int s = line.find("SN");
			int e = line.find("LN");
			string name = "";
			for (int i = s+3; i < line.size(); i++) {
				if (line[i] == ' ' || line[i] == '\t')
					break;
				name += line[i];
			}
			string ln = "";
			for (int i = e+3; i < line.size(); i++) {
				if (line[i] == ' ' || line[i] == '\t')
					break;
				ln += line[i];
			}
			length = stoi(ln);
			sum += length;
			chr_len.insert({name, length});
			// cerr << name << " | " << length << endl;
		}
		getline(fin, line);
	}
	if (chr_len.size() == 0) {
		show_progress = false;
		Logger::instance().info("No sam header detected, progress will not be shown. Please wait until partitioning finishes.\n");
	}

	string comment = "Partitioning";
	if (show_progress)
		progress.update((0.0/(float)sum) * 100, comment);

    sam_record record(line);
    curr_chr = record.rname;
    genome_offset = record.pos - array_offset;

	uint64_t processed = 0;
    while(!fin.eof()) {
		if (show_progress)
			progress.update(((float)(processed + genome_offset)/(float)sum) * 100, comment);
        process_full = false;
        //cerr << "-----" << endl;
        auto next = fill(locs, genome_offset, fin);
        //cerr << genome_offset << " | " << array_offset << endl;
        if (next.first != curr_chr || (next.first == curr_chr && genome_offset + window < next.second)) {
            next_genome_offset = next.second - 5;
            process_full = true;
            next_chr = next.first;
        }
        array_offset = clean_up(locs, genome_offset, process_full, next.second);
        if (process_full) {
            genome_offset = next_genome_offset;
            if (curr_chr != next_chr) {
                // cerr << "processed: " << curr_chr << endl;
				processed += chr_len[curr_chr];
            }
            curr_chr = next_chr;
        }
        else
            genome_offset += window - array_offset;
    //    cerr << genome_offset << " | " << array_offset << endl;
    }

	if (show_progress)
		progress.update(100.0, comment);

//    //cerr << "-----" << endl;
//    for (int i = 0; i < locs.size(); i++) {
//        //cerr << genome_offset + i << ": " << locs[i].support << " | " << locs[i].window_support << " | " << locs[i].nbr_stretch << endl;
//    }

//    string curr_chr = "";
//    for (string line; getline(fin, line);) {
//        sam_record record(line);
//
//        if ((!record.is_fully_mapped) && ((record.flag & 2048) != 2048) && (record.mapq > MIN_Q)) {
//            if (record.rname != curr_chr) {
//                curr_chr = record.rname;
//                clean_up(locs, offset);
//                for (int i = 0; i < margin; i++) {
//                    BreakpointCandidate tmp;
//                    locs[i] = tmp;
//                }
//                offset = record.pos - margin;
//            }
//
//            if (record.pos >= offset + window) {
//                clean_up(locs, offset);
//                offset = record.pos - margin;
////                //cerr << "offset: " << offset << endl;
//            }
//
//            if (record.head_clip_range != 0) {
//                locs[record.pos - offset].support++;
//                locs[record.pos - offset].reads.push_back({record.qname, record.seq});
//                locs[record.pos - offset].chr = record.rname;
//            }
//
//            if (record.tail_clip_range != 0) {
//                locs[record.end_pos - offset].support++;
//                locs[record.end_pos - offset].reads.push_back({record.qname, record.seq});
//                locs[record.end_pos - offset].chr = record.rname;
//            }
//
//            if (record.insertion_pos != record.pos) {
//                locs[record.insertion_pos - offset].support++;
//                locs[record.insertion_pos - offset].reads.push_back({record.qname, record.seq});
//                locs[record.insertion_pos - offset].chr = record.rname;
//            }
//        }
//    }
//    clean_up(locs, offset);
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