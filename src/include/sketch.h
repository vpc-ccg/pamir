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
	int seq_id;
	int offset;
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
		std::string file;
        std::vector<std::pair<uint64_t, Location> > ref_minimizers_vec;
        std::vector<std::pair<uint64_t, int> > query_minimizers_vec;
        void update_ref_sketch();
        void update_query_sketch();

	public:
		std::map<int, std::string> sequences_names;
		std::map<uint64_t, std::vector<Location> > ref_minimizers;
		std::unordered_set<uint64_t> query_minimizers;
		
		Sketch(std::string path, int k = 15, int w = 10);
		void build_sketch();
		void build_sketch(std::vector<std::string>);
        void sketch_query(std::vector<std::string> reads, int k = 15, int w = 10);
		void get_ref_minimizers(char* read, int read_id, int len);
        void get_query_minimizers(char* read, int read_id, int len);
        std::vector<std::pair<std::string, std::pair<std::pair<int, int>, int> > > find_cuts();
        std::vector<std::pair<std::string, std::pair<std::pair<int, int>, int> > > classify_reads
                    (std::map<int, std::vector<std::pair<int, uint64_t> > > hits,std::vector<std::pair<int, cut> > cuts_tmp);
};


#endif
