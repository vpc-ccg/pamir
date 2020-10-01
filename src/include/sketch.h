#ifndef SKETCH_H
#define SKETCH_H

#include <map>
#include <string>
#include <vector>
#include <utility>
#include <cstdint>

struct Location {
	int seq_id;
	int offset;
};

class Sketch {
	private:
		int kmer_size;
		int window_size;
		std::string file;

	public:
		std::map<int, std::string> sequences;
		std::map<uint64_t, std::vector<Location> > minimizers;
		
		Sketch(std::string path, int k = 15, int w = 10);
		Sketch(std::vector<std::string> reads, int k = 15, int w = 10);
		void build_sketch();
		void build_sketch(std::vector<std::string>);
		void get_minimizers(char* read, int read_id, int len);
		void print();
		std::vector<std::pair<std::string, std::pair<int, int> > > find_cuts(Sketch q);
};


#endif
