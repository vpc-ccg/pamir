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

#include "aligner.h"
#include "progressbar.h"

const int BUFFSIZE = 2000000000;

typedef uint64_t hash_t;
typedef uint64_t mem_offset_t; //CHANGE
typedef uint32_t id_t;
typedef uint32_t hash_size_t;
typedef uint16_t offset_t;

enum orientation_en : u_int8_t { FRW = 0, REV = 1};
enum type_en : u_int8_t {BIMODAL = 0, LEFT = 1, RIGHT = 2, MISC = 3, SINGLE_PEAK = 4};

//CHECK OCCURNCES
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
        return seq_id < y.seq_id;
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

//TO DO
struct seq_data {
    string name;
    offset_t size;
};

class Sketch {
	private:
        int kmer_size = 15;
        gzFile gz_fin;
        char *zbuffer;
		int window_size = 10;
        int read_id = 0;
        int freq_th = 0;
        uint64_t file_size;
        string lr_path;
		string dat_path;
        int32_t buff_pos = 0;
        int32_t buff_size = 0;
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
		Sketch(string dat_path, int k = 15, int w = 10);
        Sketch(string lr_path, string dat_path, int k = 15, int w = 10);
        void build_sketch();
		void build_sketch_mt(int, const ProgressBar, vector<pair<hash_t, Location> > &);
        void get_ref_minimizers(char*, id_t, int, vector<pair<hash_t, Location> > &);

        vector<cut> query(vector<string>&, bool);
        void build_query_sketch(vector<string>&, vector<pair<uint64_t, int> >&, unordered_set<hash_t> &, unordered_set<hash_t> &);
        void get_query_minimizers(char*, id_t, offset_t, vector<pair<uint64_t, int> > &);

        cut find_range(vector<hit>& hits, mem_offset_t start, hash_size_t size);
        vector<cut> find_cuts(bool, unordered_set<hash_t> &, unordered_set<hash_t> &);
        void classify_reads(vector<hit> &, vector<cut> &);

        float compare_sequences(string& seq_a, string& seq_b);
        float minimizer_similarity(unordered_set<hash_t>& ref, vector<hit>& q, int start, int size);
        pair<float, float> minimizer_similarity_new(unordered_set<hash_t>& ref, vector<hit>& q, int start, int size);
        float minimizer_similarity(unordered_set<hash_t>& ref, unordered_set<hash_t>& q);
};

#endif
