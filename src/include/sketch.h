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

const int BIMODAL = 0;
const int SINGLE_PEAK = 4;
const int LEFT = 1;
const int RIGHT = 2;
const int MISC = 3;
const int FRW = 0;
const int REV = 1;

const int BUFFSIZE = 2000000000;

struct Location {
    uint32_t seq_id;
    uint16_t offset;
};

struct cut {
    pair<int, int> range;
    int type;
    pair<int, int> peak1;
    pair<int, int> peak2;
    int number_of_minimizers;
    int orientation;
};

struct hash_t {
    uint64_t hash_value;
    uint64_t offset;
    bool operator() (const uint64_t& i, const hash_t &b) {
        return (i < b.hash_value);
    }
    bool operator() (const hash_t &b, const uint64_t& i) {
        return (b.hash_value < i);
    }
};

class Sketch {
	private:
        int kmer_size;
        gzFile gz_fin;
        char *zbuffer;
		int window_size;
        int read_id = 0;
        int freq_th = 0;
        uint64_t file_size;
        string lr_path;
		string dat_path;
        int32_t buff_pos = 0;
        int32_t buff_size = 0;
		static const int thread_cnt = 16;

        vector<hash_t> hashes;
        vector<Location> ref_minimizers;
        vector<pair<string, uint16_t> > names;

        void load();
        void compute_freq_th();
        void update_query_sketch(int ort);
        void dump(vector<pair<uint64_t, Location> > &ref_minimizers_vec);

        void read_buffer();
        uint32_t read_line(string& seq);
        pair<uint64_t, uint64_t> find_hit(const uint64_t &hv);
        void merge(vector<pair<uint64_t, Location> > &a, vector<pair<uint64_t, Location> > &b);

	public:
        Sketch();
		Sketch(string dat_path, int k = 15, int w = 10);
        Sketch(string lr_path, string dat_path, int k = 15, int w = 10);
        void build_sketch();
		void build_sketch_mt(int, const ProgressBar, vector<pair<uint64_t, Location> > &);
        void get_ref_minimizers(char*, int, int, vector<pair<uint64_t, Location> > &);

        vector<pair<string, pair<pair<int, int>, pair<int, int> > > > query(vector<string>&, bool);
        void build_query_sketch(vector<string>&, vector<pair<uint64_t, int> >&, unordered_set<uint64_t> &, unordered_set<uint64_t> &);
        void get_query_minimizers(char*, int, int, vector<pair<uint64_t, int> > &);

        vector<pair<string, pair<pair<int, int>, pair<int, int> > > > find_cuts(bool, unordered_set<uint64_t>, unordered_set<uint64_t>);
        vector<pair<string, pair<pair<int, int>, pair<int, int > > > > classify_reads (map<int, vector<pair<int, uint64_t> > >, vector<pair<int, cut> >);
};

void fix_reverse(vector<pair<pair<string, string>, pair<pair<int, int>, pair<int, int> > > > &cuts);

#endif
