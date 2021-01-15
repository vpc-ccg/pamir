#ifndef __CUT_RANGES__
#define __CUT_RANGES__

#include <map>
#include <vector>

typedef std::pair<int, std::pair<int, std::string> > read_cut;

class cut_ranges {
	private:
        std::string lr_path;
		std::map<std::string, std::uint64_t> read_offsets;
		std::map<std::string, std::pair<int, int> > lr_ranges;
		std::vector<std::string> names;
        std::vector<std::pair<int, std::pair<int, int> > > ranges;

	public:
		cut_ranges();
        cut_ranges(const std::string &, bool build_index = false);
		cut_ranges(const std::string &, const std::string &, const std::string &, bool build_index = false);
        std::map<std::string, std::pair<int, std::string> > reads;
        std::vector<read_cut> reads_new;
		void add_range(int id, int start, int end);
		void extract();
        void extract_new();
		std::string get_cut(int id, int start, int end);
		std::pair<std::string, std::string> get_cut_new(int id, int start, int end);
        read_cut find_read(int id);
};

std::pair<std::string, std::pair<int, int>> cut_consensus_bimodal(std::vector<std::string> alignments,
                                                                  int left_reads, int right_reads, int bimodal_reads);
std::pair<std::string, std::pair<int, int>> cut_consensus_single(std::vector<std::string> alignments);
std::vector<std::vector<std::string> > cluster_reads(std::vector<std::string> msa, int l, int r);

#endif