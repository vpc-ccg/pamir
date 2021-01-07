#ifndef SKETCH_H
#define SKETCH_H

#include <map>
#include <string>
#include <vector>
#include <utility>
#include <cstdint>
#include <unordered_set>

const int BIMODAL = 0;
const int SINGLE_PEAK = 4;
const int LEFT = 1;
const int RIGHT = 2;
const int MISC = 3;

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
};

class Sketch {
	private:
		int kmer_size;
		int window_size;
        std::string lr_path;
		std::string dat_path;
        std::vector<std::pair<uint64_t, Location> > ref_minimizers_vec;
        std::vector<std::pair<uint64_t, int> > query_minimizers_vec;
        void update_ref_sketch();
        void update_query_sketch();
        void dump();
        void load();

	public:
        std::map<uint32_t , std::pair<std::string, uint16_t> > sequences;
		std::map<uint64_t, std::vector<Location> > ref_minimizers;
		std::unordered_set<uint64_t> query_minimizers;

		Sketch();
		Sketch(std::string dp);
        Sketch(std::string lp, std::string dp, int k = 15, int w = 10);
		void build_sketch();
		void build_sketch(std::vector<std::string>);
        void sketch_query(std::vector<std::string> reads, int k = 15, int w = 10);
        void sketch_query(std::string query, int k = 15, int w = 10);
		void get_ref_minimizers(char* read, int read_id, int len);
        void get_query_minimizers(char* read, int read_id, int len);
        std::vector<std::pair<std::string, std::pair<std::pair<int, int>, int> > > find_cuts(bool classify);
        std::vector<std::pair<std::string, std::pair<std::pair<int, int>, int> > > classify_reads
                    (std::map<int, std::vector<std::pair<int, uint64_t> > > hits,std::vector<std::pair<int, cut> > cuts_tmp);
};


#endif
