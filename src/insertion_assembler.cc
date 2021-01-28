#include "insertion_assembler.h"

#include <mutex>
#include <string>
#include <iterator>
#include <algorithm>
#include "common.h"

using namespace std;

mutex extract_lock;

InsertionAssembler::InsertionAssembler(Sketch* sketch, cut_ranges* read_extractor) {
    lr_sketch = sketch;
    extractor = read_extractor;
}

string InsertionAssembler::build_segment(vector<string> &cuts) {
    spoa::Graph graph{};
    graph.Clear();
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 2, -32, -64, -1);
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

vector<string> InsertionAssembler::extract_reads(map<id_t, cut> &cuts) {
    vector<string> ext;
    for (auto it = cuts.begin(); it != cuts.end(); it++) {
        extract_lock.lock();
        cerr << "------" << endl;
        cerr << it->first << " | " << it->second.range.start << " - " << it->second.range.end << endl;
        string r = extractor->get_cut(lr_sketch->sequences[it->first].first, it->second.range.start, it->second.range.end);
        extract_lock.unlock();
        if (it->second.orientation == REV)
            ext.push_back(reverse_complement(r));
        else
            ext.push_back(r);
    }
    sort(ext.begin(), ext.end(), [](const auto &a, const auto &b) { return a.size() > b.size();});
    return ext;
}

map<id_t, cut> check_end(map<id_t, cut> &l, map<id_t, cut> &r) {
    map<id_t, cut> mid;
    for (auto it = l.begin(); it != l.end(); it++) {
        auto f = r.find(it->first);
        if (f != r.end() && it->second.orientation == f->second.orientation) {
            offset_t start, end;
            id_t id = it->first;
            if (it->second.orientation == FRW) {
                start = f->second.range.start;
                end = it->second.range.end;
            }
            else {
                start = it->second.range.start;
                end = f->second.range.end;
            }
            range_s range = {start, end};
            cut merged_cut;
            merged_cut.range = range;
            merged_cut.seq_id = id;
            merged_cut.orientation = FRW;
            mid.insert({id, merged_cut});
        }
    }
    return mid;
}

//TODO: Remove map
map<id_t, cut> InsertionAssembler::find_cuts(string& segment, bool left) {
    vector<string> dummy;
    dummy.push_back(segment);
    vector<cut> cuts = lr_sketch->query(dummy, false);
    map<id_t, cut> ans;
    for (auto it = cuts.begin(); it != cuts.end(); it++) {
        offset_t start, end;
        if (it->orientation == FRW) {
            start = left ? it->range.start : 1;
            end = left ? lr_sketch->sequences[it->seq_id].second - 1: it->range.end;
        }
        else {
            start = left ? 1 : it->range.start;
            end = left ? it->range.end : lr_sketch->sequences[it->seq_id].second - 1;
        }
        it->range.start = start;
        it->range.end = end;
        ans.insert({it->seq_id, *it});
    }

    return ans;
}

string InsertionAssembler::get_overlap(vector<string>& l, vector<string>& r, string& m) {
    spoa::Graph graph{};
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

pair<string, int> InsertionAssembler::assemble(vector<string>& left_reads, vector<string>& right_reads, string &output) {
	output += "--- Building Long Insertion ---";

    string left_seg, right_seg, lanchor, ranchor;
    vector<string> lsegs, rsegs;
    map<id_t, cut> lcuts, rcuts, mid;
    unordered_set<hash_t> l_minimizers_frw, l_minimizers_rev, r_minimizers_frw, r_minimizers_rev;
    int cut_size, support = left_reads.size() + right_reads.size();
    int step = 1;

    string l_tmp, r_tmp;

    while (true) {
        left_seg = build_segment(left_reads);
        lsegs.push_back(left_seg);
        cut_size = min((int) left_seg.size() / 2, 300);
        if (step != 0)
            lanchor = left_seg.substr(left_seg.size() - cut_size, cut_size);
        else
            lanchor = left_seg;
		lcuts = find_cuts(lanchor, true);

        right_seg = build_segment(right_reads);
        rsegs.push_back(right_seg);
        cut_size = min((int) right_seg.size() / 2, 300);
        if (step != 0)
            ranchor = right_seg.substr(0, cut_size);
        else
            ranchor = right_seg;
        rcuts = find_cuts(ranchor, false);

        mid = check_end(lcuts, rcuts);
        if (mid.size() != 0) {
			output += "\n===> Insertion built in "+ to_string(step) + " steps.";
			break;
		}
		else if (step == 10) {
			output += "\n===> Could not build insertion properly in " + to_string(step) + " steps.";
			break;
		}

        support += left_reads.size() + right_reads.size();
        left_reads = extract_reads(lcuts);
        right_reads = extract_reads(rcuts);

		output += "\n* " + to_string(step) + ":\n";
		output += "          Left Segment :  " + left_seg + "\n";
		output += "          Right Segment:  " + right_seg + "\n";
		output += "          Picked " + to_string(left_reads.size()) + " new left cuts\n";
		output += "          Picked " + to_string(right_reads.size()) + " new right cuts\n";

        step += 1;
    }

    vector<string> middle_reads = extract_reads(mid);
    string middle = build_segment(middle_reads);
    support += middle_reads.size();

    string insertion = get_overlap(lsegs, rsegs, middle);

    return {insertion, support};
}