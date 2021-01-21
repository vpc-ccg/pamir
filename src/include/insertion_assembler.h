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

        string build_segment(vector<string>& cuts);
        vector<string> extract_reads(map<id_t, cut>& cuts);
        map<id_t, cut> find_cuts(string&, bool left);
        string get_overlap(vector<string>& l, vector<string>& r, string& m);
    public:
        InsertionAssembler(Sketch* sketch, cut_ranges* read_extractor);
        pair<string, int> assemble(vector<string>& left_reads, vector<string>& right_reads, string &output);
};


#endif
