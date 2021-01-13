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

vector<string> InsertionAssembler::extract_reads(map<string, pair<pair<int, int>, int> > &cuts) {
    vector<string> ext;
    for (auto it = cuts.begin(); it != cuts.end(); it++) {
        extract_lock.lock();
        string r = extractor->get_cut(it->first, it->second.first.first, it->second.first.second);
        extract_lock.unlock();
        if (it->second.second == REV)
            ext.push_back(reverse_complement(r));
        else
            ext.push_back(r);
    }
    sort(ext.begin(), ext.end(), [](const auto &a, const auto &b) { return a.size() > b.size();});
    return ext;
}

map<string, pair<pair<int, int>, int> > check_end(map<string, pair<pair<int, int>, int> > &l, map<string, pair<pair<int, int>, int> > &r) {
    map<string, pair<pair<int, int>, int> > mid;
    for (auto it = l.begin(); it != l.end(); it++) {
        auto f = r.find(it->first);
        if (f != r.end() && it->second.second == f->second.second) {
            if (it->second.second == FRW)
                mid.insert({it->first, {{f->second.first.first, it->second.first.second}, it->second.second}});
            else
                mid.insert({it->first, {{it->second.first.first, f->second.first.second}, it->second.second}});
        }
    }
    return mid;
}

map<string, pair<pair<int, int>, int> > InsertionAssembler::find_cuts(string& segment, bool left) {
    vector<string> dummy;
    dummy.push_back(segment);
    vector<pair<string, pair<pair<int, int>, pair<int, int> > > > cuts = lr_sketch->query(dummy, false);
    map<string, pair<pair<int, int>, int> > ans;
    for (auto it = cuts.begin(); it != cuts.end(); it++) {
        int start, end;
        if (it->second.second.second == FRW) {
            start = left ? it->second.first.first : 1;
            end = left ? it->second.second.first : it->second.first.second;
        }
        else {
            start = left ? 1 : it->second.first.first;
            end = left ? it->second.first.second : it->second.second.first;
        }
        ans.insert({it->first, {{start, end}, it->second.second.second}});
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

pair<string, int> InsertionAssembler::assemble(vector<string>& left_reads, vector<string>& right_reads) {
	Logger::instance().info("--- Building Long Insertion ---");

    string left_seg, right_seg, lanchor, ranchor;
    vector<string> lsegs, rsegs;
    map<string, pair<pair<int, int>, int> > lcuts, rcuts, mid;
    unordered_set<uint64_t> l_minimizers_frw, l_minimizers_rev, r_minimizers_frw, r_minimizers_rev;
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
			Logger::instance().info("\n===> Insertion built in %d steps.", step);
			break;
		}
		else if (step == 10) {
			Logger::instance().info("\n===> Could not build insertion properly in %d steps.", step);
			break;
		}

        support += left_reads.size() + right_reads.size();
        left_reads = extract_reads(lcuts);
        right_reads = extract_reads(rcuts);

		Logger::instance().info("\n* STEP %d:\n", step);
		Logger::instance().info("          Left Segment :  %s\n", left_seg.c_str());
		Logger::instance().info("          Right Segment:  %s\n", right_seg.c_str());
		Logger::instance().info("          Picked %d new left cuts\n", left_reads.size());
		Logger::instance().info("          Picked %d new right cuts\n", right_reads.size());

        step += 1;
    }

    vector<string> middle_reads = extract_reads(mid);
    string middle = build_segment(middle_reads);
    support += middle_reads.size();

    string insertion = get_overlap(lsegs, rsegs, middle);

    return {insertion, support};
}