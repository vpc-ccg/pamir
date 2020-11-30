#include "insertion_assembler.h"

#include <string>
#include <cstring>
#include <iterator>
#include <algorithm>

using namespace std;

InsertionAssembler::InsertionAssembler(const string &dat_path, const string &lr_path) : lr_path(lr_path) {
    lr_sketch = Sketch(dat_path, true);
}

string InsertionAssembler::build_segment(vector<string> cuts) {
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 2, -32, -64, -1);
    graph.Clear();
    for (int i = 0; i < cuts.size(); i++) {
        auto alignment = alignment_engine->Align(cuts[i], graph);
        graph.AddAlignment(alignment, cuts[i]);
    }
    vector<string> msa = graph.GenerateMultipleSequenceAlignment(true);

    int th = 0.7 * (msa.size() - 1);

    string consensus = msa[msa.size() - 1];

    int left = 0, right = 0;
    int gap_cnt;
    for (int idx = 0; idx < msa[0].size(); idx = idx + 50) {
        gap_cnt = 0;
        for (int i = 0; i < msa.size() - 1; i++) {
            if (msa[i][idx] == '-')
                gap_cnt++;
        }
        if (gap_cnt < th) {
            left = idx;
            break;
        }
    }

    for (int idx = 0; idx < msa[0].size(); idx = idx+50) {
        gap_cnt = 0;
        for (int i = 0; i < msa.size() - 1; i++) {
            if (msa[i][msa[i].size() - idx - 1] == '-')
                gap_cnt++;
        }
        if (gap_cnt < th) {
            right = idx;
            break;
        }
    }

    string cut = consensus.substr(left, consensus.size() - right - left);
    cut.erase(std::remove(cut.begin(), cut.end(), '-'), cut.end());
    return cut;
}

vector<string> InsertionAssembler::extract_reads(map<string, pair<int, int> > cuts) {
    cut_ranges read_extractor = cut_ranges(lr_path);
    for (auto it = cuts.begin(); it != cuts.end(); it++) {
        read_extractor.add_range(it->first, it->second.first, it->second.second);
    }
    read_extractor.extract();
    vector<string> ext;
    for (auto it = cuts.begin(); it != cuts.end(); it++) {
        ext.push_back(read_extractor.get_cut(it->first, it->second.first, it->second.second));
    }
    sort(ext.begin(), ext.end(), [](const auto &a, const auto &b) { return a.size() > b.size();});
    return ext;
}

map<string, pair<int, int> > check_end(map<string, pair<int, int> > l, map<string, pair<int, int> > r) {
    map<string, pair<int, int> > mid;
    for (auto it = l.begin(); it != l.end(); it++) {
        auto f = r.find(it->first);
        if (f != r.end()) {
            mid.insert({it->first, {f->second.first, it->second.second}});
        }
    }
    return mid;
}

map<string, pair<int, int> > InsertionAssembler::find_cuts(bool left) {
    vector<pair<string, pair<pair<int, int>, int> > > cuts = lr_sketch.find_cuts(false);
    map<string, pair<int, int> > ans;
    for (auto it = cuts.begin(); it != cuts.end(); it++) {
        if (!(it->second.first.first > 300 && it->second.first.second < it->second.second - 300 ))
            continue;
        int start = left ? it->second.first.first : 1;
        int end = left ? it->second.second : it->second.first.second;
        ans.insert({it->first, {start, end}});
    }
    return ans;
}

string InsertionAssembler::get_overlap(vector<string> l, vector<string> r, string m) {
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kOV, 2, -32, -64, -1);
    graph.Clear();
    for (int i = 0; i < l.size(); i++) {
        auto alignment = alignment_engine->Align(l[i], graph);
        graph.AddAlignment(alignment, l[i]);
    }
    auto alignment = alignment_engine->Align(m, graph);
    graph.AddAlignment(alignment, m);
    for (int i = r.size() - 1; i >= 0; i--) {
        auto alignment = alignment_engine->Align(r[i], graph);
        graph.AddAlignment(alignment, r[i]);
    }

    auto msa = graph.GenerateMultipleSequenceAlignment(true);
    string consensus = msa[msa.size() - 1];
    consensus.erase(std::remove(consensus.begin(), consensus.end(), '-'), consensus.end());

    return consensus;
}

pair<string, int> InsertionAssembler::assemble(vector<string> left_reads, vector<string> right_reads) {
    string left_seg, right_seg, lanchor, ranchor;
    vector<string> lsegs, rsegs;
    map<string, pair<int, int> > lcuts, rcuts, mid;
    int cut_size, support = left_reads.size() + right_reads.size();

    while (true) {
        left_seg = build_segment(left_reads);
        lsegs.push_back(left_seg);
        cut_size = min((int) left_seg.size() / 2, 300);
        lanchor = left_seg.substr(left_seg.size() - cut_size, cut_size);
        lr_sketch.sketch_query(lanchor);
        lcuts = find_cuts(true);

        right_seg = build_segment(right_reads);
        rsegs.push_back(right_seg);
        cut_size = min((int) right_seg.size() / 2, 300);
        ranchor = right_seg.substr(0, cut_size);
        lr_sketch.sketch_query(ranchor);
        rcuts = find_cuts(false);

        mid = check_end(lcuts, rcuts);
        if (mid.size() != 0)
            break;

        support += left_reads.size() + right_reads.size();
        left_reads = extract_reads(lcuts);
        right_reads = extract_reads(rcuts);
    }

    vector<string> middle_reads = extract_reads(mid);
    string middle = build_segment(middle_reads);
    support += middle_reads.size();

    string insertion = get_overlap(lsegs, rsegs, middle);

    return {insertion, support};
}