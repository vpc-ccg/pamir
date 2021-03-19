#ifndef SKETCH_H
#define SKETCH_H

#include <map>
#include <unordered_map>
#include <zlib.h>
#include <string>
#include <vector>
#include <utility>
#include <cstdint>
#include <fstream>
#include <unordered_set>

#include "chain.h"
#include "aligner.h"
#include "progressbar.h"

const int BUFFSIZE = 2000000000;

const int GENOME_ANCHOR_LEN = 100;
const int GENOME_ANCHOR_CUTOFF = 50;
//TODO: = k + 2*w
const int BOUNDARY_DISTANCE_CUTOFF = 100;
const float INSERTION_MINIMIZERS_CUTOFF = 0.125;
const float TCK = 0.9;

typedef uint32_t id_t;
typedef uint64_t hash_t;
typedef uint16_t offset_t;
typedef uint32_t hash_size_t;
typedef uint64_t mem_offset_t;
typedef pair<hash_t, int> minimizer_t;

enum orientation_en : u_int8_t {FRW = 0, REV = 1};
enum type_en : u_int8_t {BIMODAL = 0, PARTIAL_LEFT = 1, PARTIAL_RIGHT = 2, OVERLAPPING_READ = 3, DROPPED = 4};

struct cut_stats {
    int bimodal_cnt = 0;
    int bimodal_sum = 0;
    int overlapping_cnt = 0;
    vector<int> left_cuts_size;
    vector<int> right_cuts_size;
};

struct Location {
    id_t seq_id;
    offset_t offset;

    bool operator <(const Location &y) {
        if (seq_id < y.seq_id)
            return true;
        else if (seq_id == y.seq_id) {
            return offset < y.offset;
        }
        return false;
    }
};

struct range_s {
    offset_t start;
    offset_t end;
};

struct hit {
    id_t seq_id;
    hash_t hash_value;
    mem_offset_t offset;
    mem_offset_t genome_offset;

    bool operator <(const hit &y) {
        return seq_id == y.seq_id ? hash_value < y.hash_value : seq_id < y.seq_id;
    }

    bool operator () (const id_t &i, const hit &a) const {
        return i < a.seq_id;
    }

    bool operator () (const hit &a, const id_t &i) const {
        return a.seq_id < i;
    }

    bool operator() (const hit &a, const hit &b) {
        return a.offset < b.offset;
    }
};

struct cut {
    id_t seq_id;
    type_en type;
    range_s range;
    hash_size_t size;
    range_s genome_range;
    int breakpoint_distance;
    orientation_en orientation;
    int estimated_insertion = -1;

    bool operator <(const cut &b) {
        if (type < b.type) {
            return true;
        }
        else if (type == b.type) {
            if (range.end - range.start > b.range.end - b.range.start)
                return true;
            return false;
        }
        else
            return false;
    }
};

struct minimizer {
    hash_t hash_value;
    mem_offset_t offset;

    bool operator() (const minimizer &b, const hash_t &i) {
        return (b.hash_value < i);
    }
};

class Sketch {
	private:

        gzFile gz_fin;
        char *zbuffer;
        string lr_path;
        int freq_th = 0;
        int read_id = 0;
		string dat_path;
        int kmer_size = 15;
        uint64_t file_size;
		int window_size = 10;
        int32_t buff_pos = 0;
        int32_t buff_size = 0;
		static const int thread_cnt = 4;

        vector<minimizer> minimizers;
        vector<Location> ref_minimizers;

        void load();
        void compute_freq_th();
        void dump(vector<pair<hash_t, Location> > &ref_minimizers_vec);

        void read_buffer();
        uint32_t read_line(string& seq);
        pair<mem_offset_t, hash_size_t> find_hit(const hash_t &hv);
        void merge(vector<pair<hash_t, Location> > &a, vector<pair<hash_t, Location> > &b);

	public:
        vector<pair<string, offset_t> > sequences;
        ClaspChain claspChain;

        Sketch();
		Sketch(string dat_path, int len, int distance, int k = 15, int w = 10);
        Sketch(string lr_path, string dat_path, int k = 15, int w = 10);
        void build_sketch();
		void build_sketch_mt(int, const ProgressBar, vector<pair<hash_t, Location> > &);
        void get_ref_minimizers(char*, id_t, int, vector<pair<hash_t, Location> > &);
        pair<vector<cut>, int> query(vector<string>&, bool, string&, string&);
        void get_query_minimizers(char*, id_t, offset_t, vector<pair<uint64_t, int> > &);
        pair<vector<cut>, int> find_cuts_all(string& ref_l, string& ref_r, bool long_insertion, vector<string>& reads);
        void get_genome_hits_new(string& ref, vector<minimizer_t>& minimizers, vector<hit>& candidates);
        vector<hit> get_hits_new(vector<pair<uint64_t, int> >& query_frw);
        vector<cut> find_cuts_with_chain(string ref_l, string ref_r, cut_stats&, orientation_en, unordered_set<id_t>&, bool);
        vector<seed> create_seeds_new(vector<hit>& hits, int start, int size);
        MaxChainInfo get_genome_anchor_new_new(vector<hit>& left_anchor_hits, int start, int size);
        void find_left_cuts(vector<hit>& read_candidates, vector<cut>&, int, orientation_en, int genome_minimizers_cnt, unordered_set<id_t>& insertion_minimizers);
        void find_right_cuts(vector<hit>& read_candidates, vector<cut>&, int, orientation_en, int genome_minimizers_cnt, unordered_set<id_t>& insertion_minimizers);
        void merge_candidates(vector<cut>& left_candidates, vector<cut>& right_candidates, vector<cut>& merged_candidates,
                          cut_stats&, bool);
        void get_insertion_minimizers_new(vector<string>& short_reads, string& genome_left, string& genome_right, unordered_set<id_t>& candidates);
        void get_unique_minimizers(vector<string> &reads, unordered_set<hash_t>& insertion_minimizers);
        inline int max_kmer_count(MaxChainInfo chain, int anchor_length, int read_length, int);

        double get_minimizers_time = 0;
        double find_cuts_with_chain_time = 0;
        double insertion_estimation_time = 0;
        double readjustment_time = 0;
        double get_minimizers_time_p1 = 0;
        double get_minimizers_time_p2 = 0;
        double get_minimizers_time_p3 = 0;
        double get_minimizers_time_p4 = 0;
        double get_minimizers_time_p0 = 0;

        double hits_time = 0;
        double finding_time = 0;
        double merging_time = 0;
        double chaining_time = 0;
        double clasp_time = 0;
        double seed_time = 0;
        int tmp_cnt = 0;
};

#endif
