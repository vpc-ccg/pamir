#ifndef SKETCH_H
#define SKETCH_H

#include <map>
#include <string>
#include <vector>
#include <utility>
#include <cstdint>
#include <unordered_set>

struct Location {
	int seq_id;
	int offset;
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
        std::vector<std::pair<std::string, std::pair<int, int> > > find_cuts();
};


#endif
