#ifndef __CUT_RANGES__
#define __CUT_RANGES__

#include <map>
#include <vector>

#include "sketch.h"

typedef std::pair<id_t, std::pair<offset_t, std::string> > read_cut;

class cut_ranges {
	private:
        std::string lr_path;
        ifstream lr_file;
		std::map<std::string, uint64_t> read_offsets;
		std::vector<std::string> names;
        std::vector<std::pair<id_t, range_s> > ranges;

	public:
		cut_ranges();
        cut_ranges(const std::string &, bool build_index = false);
		cut_ranges(const std::string &, const std::string &, const std::string &, bool build_index = false);
        std::vector<read_cut> reads;
		void add_range(id_t id, offset_t start, offset_t end);
        void extract();
        std::string get_cut(std::string id, offset_t start, offset_t end);
		std::pair<std::string, std::string> get_cut(id_t id, offset_t start, offset_t end);
        read_cut find_read(id_t id);
};

std::pair<std::string, std::pair<int, int>> cut_consensus_bimodal(std::vector<std::string> alignments,
                                                                  int left_reads, int right_reads, int bimodal_reads);
std::pair<std::string, std::pair<int, int>> cut_consensus_single(std::vector<std::string> alignments);
std::vector<std::vector<std::string> > cluster_reads(std::vector<std::string> msa, int l, int r);

#endif