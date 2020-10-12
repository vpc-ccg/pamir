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

std::string cut_consensus(std::vector<std::string> alignments);

#endif