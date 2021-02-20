#ifndef SKETCH_H
#define SKETCH_H

#include <map>
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

typedef uint64_t hash_t;
typedef uint64_t mem_offset_t;
typedef uint32_t id_t;
typedef uint32_t hash_size_t;
typedef uint16_t offset_t;

enum orientation_en : u_int8_t { FRW = 0, REV = 1};
enum type_en : u_int8_t {BIMODAL = 0, PARTIAL_LEFT = 1, PARTIAL_RIGHT = 2, MISC = 3, SINGLE_PEAK = 4, LONG_INSERTION = 5};
enum anchor_en : u_int8_t {FULLY_ANCHORED = 0, LEFT_ANCHORED = 1, RIGHT_ANCHORED = 2};

struct GenomeAnchor {
    int left_frw;
    int left_rev;
    int right_frw;
    int right_rev;

    bool operator <(const int &y) {
        return (left_frw + right_frw < y) && (left_rev + right_rev < y);
    }

    type_en type() {
        int max_l = max(left_frw, left_rev);
        int max_r = max(right_frw, right_rev);

        return max_l >= max_r ? PARTIAL_LEFT : PARTIAL_RIGHT;
    }
};

struct GenomeAnchorAbs {
    int left_cnt;           //Number of common minimizers between the left genomic anchor and the read
    int right_cnt;          //Number of common minimizers between the right genomic anchor and the read

    int left_genome_cnt;
    int right_genome_cnt;

    bool is_left_anchored() {
        return left_genome_cnt != 0 && left_cnt >= 0.25 * left_genome_cnt;
    }

    bool is_right_anchored() {
        return right_genome_cnt != 0 && right_cnt >= 0.25 * right_genome_cnt;
    }

    bool is_anchored() {
        return (left_genome_cnt != 0 && left_cnt >= 0.25 * left_genome_cnt || right_genome_cnt != 0 && right_cnt >= 0.25 * right_genome_cnt);
    }

    anchor_en type() {
        if (is_left_anchored() && is_right_anchored())
            return FULLY_ANCHORED;
        else if (is_left_anchored())
            return LEFT_ANCHORED;
        else
            return RIGHT_ANCHORED;
    }
};

struct Location {
    id_t seq_id;
    offset_t offset;
};

struct range_s {
    offset_t start;
    offset_t end;
};

struct hit {
    id_t seq_id;
    mem_offset_t offset;
    hash_t hash_value;

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
    range_s range;
    type_en type;
    range_s peak1;
    range_s peak2;
    hash_size_t size;
    orientation_en orientation;
    int estimated_insertion = -1;

    bool operator <(const cut &b) {
        if (type < b.type)
            return true;
        else if (type == b.type)
            return size > b.size;
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

struct seq_data {
    string name;
    offset_t size;
};

class Sketch {
	private:
        ClaspChain claspChain;
        gzFile gz_fin;
        char *zbuffer;
        string lr_path;
        int freq_th = 0;
        int read_id = 0;
		string dat_path;
        int short_read_len;
        int kmer_size = 15;
        uint64_t file_size;
		int window_size = 10;
        int32_t buff_pos = 0;
        int32_t buff_size = 0;
		int genome_anchor_distance;
		static const int thread_cnt = 16;

        vector<minimizer> minimizers;
        vector<Location> ref_minimizers;

        void load();
        void compute_freq_th();
        void update_query_sketch(int ort);
        void dump(vector<pair<hash_t, Location> > &ref_minimizers_vec);

        void read_buffer();
        uint32_t read_line(string& seq);
        pair<mem_offset_t, hash_size_t> find_hit(const hash_t &hv);
        void merge(vector<pair<hash_t, Location> > &a, vector<pair<hash_t, Location> > &b);

	public:
        vector<pair<string, offset_t> > sequences;

        Sketch();
		Sketch(string dat_path, int len, int distance, int k = 15, int w = 10);
        Sketch(string lr_path, string dat_path, int k = 15, int w = 10);
        void build_sketch();
		void build_sketch_mt(int, const ProgressBar, vector<pair<hash_t, Location> > &);
        void get_ref_minimizers(char*, id_t, int, vector<pair<hash_t, Location> > &);

        vector<cut> query(vector<string>&, bool);
        vector<cut> query(vector<string>&, bool, string&, string&);
        void build_query_sketch(vector<string>&, vector<pair<uint64_t, int> >&, unordered_set<hash_t> &, unordered_set<hash_t> &);
        void get_query_minimizers(char*, id_t, offset_t, vector<pair<uint64_t, int> > &);

        GenomeAnchorAbs get_genome_anchor(pair<vector<hit>, vector<hit> > ref_l_hits, pair<vector<hit>, vector<hit> > ref_r_hits,
                                          id_t id, orientation_en orientation, int l_cnt, int r_cnt, vector<pair<uint64_t, int> > genome_l,
                                          vector<pair<uint64_t, int> > genome_r);
        pair<vector<hit>, vector<hit> > get_hits(unordered_set<hash_t>& query_frw, unordered_set<hash_t>& query_rev);
        cut find_range(vector<hit>& hits, mem_offset_t start, hash_size_t size);
        vector<cut> find_cuts(bool, unordered_set<hash_t> &, unordered_set<hash_t> &);
        vector<cut> find_cuts(bool, unordered_set<hash_t> &, unordered_set<hash_t> &, string ref_l, string ref_r);
        void classify_reads(vector<hit> &, vector<cut> &);

        float compare_sequences(string& seq_a, string& seq_b);
        pair<float, float> minimizer_similarity(unordered_set<hash_t>& ref, vector<hit>& q, int start, int size);
        float minimizer_similarity_single(unordered_set<hash_t>& ref, unordered_set<hash_t>& q);

        pair<vector<hit>, vector<hit> > get_hits(vector<pair<uint64_t, int> >& query_frw, vector<pair<uint64_t, int> >& query_rev);
        void get_genome_hits(string& ref_l, string& ref_r, vector<pair<uint64_t, int> >&,
                             vector<pair<uint64_t, int> >&, vector<pair<uint64_t, int> >&, vector<pair<uint64_t, int> >&,
                             pair<vector<hit>, vector<hit> >& l_hits, pair<vector<hit>, vector<hit> >& r_hits);
        float estimate_insertion(vector<hit> hits_l, vector<hit> hits_r, vector<pair<uint64_t, int> > l_minimizers,
                                     vector<pair<uint64_t, int> > r_minimizers, id_t id, GenomeAnchorAbs anchor);

        void get_insertion_minimizers(unordered_set<hash_t>& reads_frw, unordered_set<hash_t>& reads_rev, string& genome,
                                  unordered_set<hash_t>& insertion_minimizers_frw, unordered_set<hash_t>& insertion_minimizers_rev);

        vector<seed> create_seeds(vector<hit> hits, id_t id, vector<pair<uint64_t, int> > genome);
};

#endif
