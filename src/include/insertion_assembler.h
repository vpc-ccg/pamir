#ifndef INSERTIONASSEMBLER_H
#define INSERTIONASSEMBLER_H

#include <memory>
#include <unordered_set>

#include "sketch.h"
#include "cut_ranges.h"
#include "spoa/spoa.hpp"

class InsertionAssembler {
    private:
        Sketch* lr_sketch;
        cut_ranges* extractor;

        std::string build_segment(std::vector<std::string>& cuts);
        std::vector<std::string> extract_reads(std::map<std::string, std::pair<std::pair<int, int>, int > >& cuts);
        std::map<std::string, std::pair<std::pair<int, int>, int> > find_cuts(string&, bool left);
        std::string get_overlap(std::vector<std::string>& l, std::vector<std::string>& r, std::string& m);
    public:
        InsertionAssembler(Sketch* sketch, cut_ranges* read_extractor);
        std::pair<std::string, int> assemble(std::vector<std::string>& left_reads, std::vector<std::string>& right_reads);
};


#endif
