#ifndef __CUT_RANGES__
#define __CUT_RANGES__

#include <map>
#include <vector>

class cut_ranges {
	private:
		std::map<std::string, std::pair<int, int> > lr_ranges;

	public:
		std::map<std::string, std::pair<int, std::string> > reads;

		void add_range(std::string lr, int start, int end);
		void extract_reads(std::string lr_path);
		std::string get_cut(std::string lr_name, int start, int end);
};

std::pair<std::string, std::pair<int, int>> cut_consensus_bimodal(std::vector<std::string> alignments,
                                                                  int left_reads, int right_reads, int bimodal_reads);
std::pair<std::string, std::pair<int, int>> cut_consensus_single(std::vector<std::string> alignments);
std::vector<std::vector<std::string> > cluster_reads(std::vector<std::string> msa, int l, int r);

#endif