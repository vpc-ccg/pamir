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
    std::pair<int, int> range;
    int type;
    std::pair<int, int> peak1;
    std::pair<int, int> peak2;
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
		int window_size;
        std::string lr_path;
		std::string dat_path;
		uint64_t file_size;
		int read_id = 0;
		static const int thread_cnt = 16;
		gzFile gz_fin;
        char *zbuffer;
        int32_t buff_pos = 0;
        int32_t buff_size = 0;

        int freq_th = 0;
        std::vector<std::pair<uint64_t, int> > query_minimizers_vec;
        void update_query_sketch(int ort);
        void compute_freq_th();
        void dump(std::vector<std::pair<uint64_t, Location> > &ref_minimizers_vec);
        void load();

        void merge(std::vector<std::pair<uint64_t, Location> > &a, std::vector<std::pair<uint64_t, Location> > &b);
        void read_buffer();
        uint32_t read_line(std::string& seq);
        std::pair<uint64_t, uint64_t> find_hit(const uint64_t &hv);

	public:
        std::vector<std::pair<std::string, uint16_t> > names;
        std::vector<hash_t> hashes;
        std::vector<Location> ref_minimizers;
		std::unordered_set<uint64_t> query_minimizers;
        std::unordered_set<uint64_t> rev_query_minimizers;

        Sketch();
		Sketch(std::string dat_path);
        Sketch(std::string lr_path, std::string dat_path, int k = 15, int w = 10);
		void build_sketch(int id, const ProgressBar progress, std::vector<std::pair<uint64_t, Location> > &ref_minimizers_vec);
		void build_sketch_mt();
		void build_sketch_query(std::vector<std::string>);
        void sketch_query(std::vector<std::string> reads, int k = 15, int w = 10);
        void sketch_query(std::string query, int k = 15, int w = 10);
		void get_ref_minimizers(char* read, int read_id, int len, std::vector<std::pair<uint64_t, Location> > &ref_minimizers_vec);
        void get_query_minimizers(char* read, int read_id, int len);
        std::vector<std::pair<std::string, std::pair<std::pair<int, int>, std::pair<int, int> > > > find_cuts(bool classify);
        std::vector<std::pair<std::string, std::pair<std::pair<int, int>, std::pair<int, int > > > > classify_reads
                    (std::map<int, std::vector<std::pair<int, uint64_t> > > hits,std::vector<std::pair<int, cut> > cuts_tmp);
};

void fix_reverse(std::vector<std::pair<std::pair<std::string, std::string>, std::pair<std::pair<int, int>, std::pair<int, int> > > > &cuts);

#endif
