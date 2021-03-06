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
//    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 2, -32, -64, -1);
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 10, -2, -15, -7);
    for (int i = 0; i < cuts.size(); i++) {
        auto alignment = alignment_engine->Align(cuts[i], graph);
        graph.AddAlignment(alignment, cuts[i]);
    }
    vector<string> msa = graph.GenerateMultipleSequenceAlignment(true);

    pair<string, pair<int, int>> cut_range = cut_consensus_single(msa, 0.5);
    int left = cut_range.second.first;
    int right = cut_range.second.second;

    string consensus = msa[msa.size() - 1];

#ifdef DEBUG
    Logger::instance().debug("- MSA\n");
    for (int i = 0; i < msa.size(); i++) {
        string tmp = msa[i].insert(left, "|");
        string tmp2 = tmp.insert(right + 1, "|");
        Logger::instance().debug("%s\n", tmp2.c_str());
    }
    Logger::instance().debug("\n");
#endif

    string cut = consensus.substr(left, right - left);
    cut.erase(std::remove(cut.begin(), cut.end(), '-'), cut.end());

    return cut;

//    int th = 0.7 * (msa.size() - 1);
//
//    string consensus = msa[msa.size() - 1];
//
//    int left = 0, right = 0;
//    int gap_cnt;
//    for (int idx = 0; idx < msa[0].size(); idx = idx + 50) {
//        gap_cnt = 0;
//        for (int i = 0; i < msa.size() - 1; i++) {
//            if (msa[i][idx] == '-')
//                gap_cnt++;
//        }
//        if (gap_cnt < th) {
//            left = idx;
//            break;
//        }
//    }
//
//    for (int idx = 0; idx < msa[0].size(); idx = idx+50) {
//        gap_cnt = 0;
//        for (int i = 0; i < msa.size() - 1; i++) {
//            if (msa[i][msa[i].size() - idx - 1] == '-')
//                gap_cnt++;
//        }
//        if (gap_cnt < th) {
//            right = idx;
//            break;
//        }
//    }

//    bool up = false, down = false;
//    int step = 200;
//    int prv = 0;
//    int left = 0, right = 0;
//    int gap_cnt;
//    for (int idx = min(200, (int)msa[0].size()/2); idx < min(200, (int)msa[0].size()/2) + 1; idx += step) {
//        gap_cnt = 0;
//        for (int i = 0; i < msa.size() - 1; i++) {
//            if (msa[i][idx] == '-')
//                gap_cnt++;
//        }
//        if (gap_cnt > th) {
//            if (down)
//                step = abs(prv - idx)/2;
//            else
//                step = idx;
//            prv = idx;
//            up = true;
//            down = false;
//        }
//        else {
//            down = true;
//            if (step != 0 && !up) {
//                step = -abs(prv - idx) / 2;
//                prv = idx;
//            } else {
//                left = idx;
//                break;
//            }
//        }
//    }
//
//    up = false;
//    step = 200;
//    prv = 0;
//    for (int idx = min(200, (int)msa[0].size()/2); idx < min(200, (int)msa[0].size()/2) + 1; idx += step) {
//        gap_cnt = 0;
//        for (int i = 0; i < msa.size() - 1; i++) {
//            if (msa[i][msa[i].size() - idx - 1] == '-')
//                gap_cnt++;
//        }
//        if (gap_cnt > th) {
//            if (down)
//                step = abs(prv - idx)/2;
//            else
//                step = idx;
//            prv = idx;
//            up = true;
//            down = false;
//        }
//        else {
//            down = true;
//            if (step != 0 && !up) {
//                step = -abs(prv - idx) / 2;
//                prv = idx;
//            } else {
//                right = idx;
//                break;
//            }
//        }
//    }

//    pair<string, pair<int, int> > cut = cut_consensus_single(msa);

//    bool up = false, down = false;
//    int left = 0, right = 0;
//    int gap_cnt;
//    int step = 200, prv = 0;
//    int curr_answer = 0;
//    ////cerr << "threshold: " << th << endl;
//    ////cerr << " ----- Left ----- " << endl;
//    for (int idx = min(200, (int)msa[0].size()/2); idx < msa[0].size(); idx += step) {
//        gap_cnt = 0;
//        for (int i = 0; i < msa.size() - 1; i++) {
//            if (msa[i][idx] == '-')
//                gap_cnt++;
//        }
//        //cerr << "......" << endl;
//        //cerr << "idx: " << idx << endl;
//        //cerr << "gap cnt: " << gap_cnt << endl;
//        if (gap_cnt > th) {
//            if (down) {
//                //cerr << "curr answer: " << curr_answer << endl;
//                if (step != 1)
//                    step = abs(curr_answer - idx)/2;
//                if (step == 0)
//                    step = 1;
//                //cerr << "down==true" << endl;
//                //cerr << "new step: " << step << endl;
//            }
//            else {
//                //cerr << "curr answer: " << curr_answer << endl;
//                if (curr_answer == 0)
//                    step = 100;
//                else if (step != 1)
//                    step = abs(curr_answer - idx)/2;
//                if (step == 0)
//                    step = 1;
////                step = idx;
//                //cerr << "down==false" << endl;
//                //cerr << "new step: " << step << endl;
//            }
//            prv = idx;
//            up = true;
//            down = false;
//            //cerr << "up true, down false" << endl;
//        }
//        else {
//            down = true;
//            //cerr << "down true" << endl;
//            if (step != 0 && step != 1) {
//                //cerr << "up==false" << endl;
//                step = -abs(prv - idx) / 2;
//                //cerr << "new step: " << step << endl;
//                prv = idx;
//                curr_answer = idx;
//            } else {
//                //cerr << "up==true" << endl;
//                //cerr << "left set to " << left << endl;
//                left = idx;
//                break;
//            }
//        }
//    }
//
//    step = 200;
//    prv = 0;
//    up = false, down = false;
//    curr_answer = 0;
//    for (int idx = min(200, (int)msa[0].size()/2); idx < msa[0].size(); idx += step) {
//        gap_cnt = 0;
//        for (int i = 0; i < msa.size() - 1; i++) {
//            if (msa[i][msa[i].size() - idx - 1] == '-')
//                gap_cnt++;
//        }
//        //cerr << "......" << endl;
//        //cerr << "idx: " << idx << endl;
//        //cerr << "gap cnt: " << gap_cnt << endl;
//        if (gap_cnt > th) {
//            if (down) {
//                //cerr << "curr answer: " << curr_answer << endl;
//                if (step != 1)
//                    step = abs(curr_answer - idx)/2;
//                if (step == 0)
//                    step = 1;
//                //cerr << "down==true" << endl;
//                //cerr << "new step: " << step << endl;
//            }
//            else {
//                //cerr << "curr answer: " << curr_answer << endl;
//                if (curr_answer == 0)
//                    step = 100;
//                else if (step != 1)
//                    step = abs(curr_answer - idx)/2;
//                if (step == 0)
//                    step = 1;
////                step = idx;
//                //cerr << "down==false" << endl;
//                //cerr << "new step: " << step << endl;
//            }
//            prv = idx;
//            up = true;
//            down = false;
//            //cerr << "up true, down false" << endl;
//        }
//        else {
//            down = true;
//            //cerr << "down true" << endl;
//            if (step != 0 && step != 1) {
//                //cerr << "up==false" << endl;
//                step = -abs(prv - idx) / 2;
//                //cerr << "new step: " << step << endl;
//                prv = idx;
//                curr_answer = idx;
//            } else {
//                //cerr << "up==true" << endl;
//                //cerr << "left set to " << left << endl;
//                right = idx;
//                break;
//            }
//        }
//    }
}

//vector<string> InsertionAssembler::extract_reads(map<id_t, cut> &cuts) {
//    vector<string> ext;
//    for (auto it = cuts.begin(); it != cuts.end(); it++) {
//        extract_lock.lock();
//        string r = extractor->get_cut(lr_sketch->sequences[it->first].first, it->second.range.start, it->second.range.end);
//        extract_lock.unlock();
////        Logger::instance().debug(">%s\n%s\n", lr_sketch->sequences[it->first].first.c_str(), r.c_str());
//        if (it->second.orientation == REV)
//            ext.push_back(reverse_complement(r));
//        else
//            ext.push_back(r);
//    }
//    sort(ext.begin(), ext.end(), [](const auto &a, const auto &b) { return a.size() > b.size();});
//    return ext;
//}

vector<string> InsertionAssembler::extract_reads(vector<cut>& cuts) {
    vector<string> ext;
    for (auto i = 0; i < cuts.size(); i++) {
        extract_lock.lock();
        string r = extractor->get_cut(lr_sketch->sequences[cuts[i].seq_id].first, cuts[i].range.start, cuts[i].range.end);
        extract_lock.unlock();

        if (cuts[i].orientation == REV)
            ext.push_back(reverse_complement(r));
        else
            ext.push_back(r);
    }
    sort(ext.begin(), ext.end(), [](const auto &a, const auto &b) { return a.size() > b.size();});
    return ext;
}

void inline get_range(int& start, int& end, int s1, int e1, int s2, int e2, int l) {
    start = max(min(s1, e2) - 100, 0);
    end = min(max(s1, e2), l);
}

map<id_t, cut> InsertionAssembler::check_end(map<id_t, cut> &l, map<id_t, cut> &r) {
    map<id_t, cut> mid;
    for (auto it = l.begin(); it != l.end(); it++) {
        auto f = r.find(it->first);
        if (f != r.end() && it->second.orientation == f->second.orientation) {
            offset_t start, end;
            id_t id = it->first;
            if (it->second.orientation == FRW) {
                start = max(min(it->second.range.start, f->second.range.end) - 100, 0);
                end = min(max(it->second.range.start, f->second.range.end) + 100, (int)lr_sketch->sequences[it->first].second);
            }
            else {
                start = max(min(it->second.range.end, f->second.range.start) - 100, 0);
                end = min(max(it->second.range.end, f->second.range.start) + 100, (int)lr_sketch->sequences[it->first].second);
//                start = it->second.range.start;
//                end = f->second.range.end;
//                start = f->second.range.start;
//                end = it->second.range.end;
//                start = max(it->second.range.end - 100, 0);
//                end = min(f->second.range.start + 100, (int)lr_sketch->sequences[it->first].second);
            }
            range_s range = {start, end};
            cut merged_cut;
            merged_cut.range = range;
            merged_cut.seq_id = id;
            merged_cut.orientation = it->second.orientation;
            mid.insert({id, merged_cut});
            Logger::instance().debug("%s: %d - %d\n", lr_sketch->sequences[it->first].first.c_str(), start, end);
        }
    }
    return mid;
}

//TODO: Remove map
//map<id_t, cut> InsertionAssembler::find_cuts(string& segment, string& anchor, bool left) {
//    vector<string> seg_vec;
//    seg_vec.push_back(segment);
//    string dummy = "";
//    vector<cut> cuts;
//    if (left)
//        cuts = lr_sketch->query(seg_vec, false, anchor, dummy);
//    else
//        cuts = lr_sketch->query(seg_vec, false, dummy, anchor);
//    map<id_t, cut> ans;
//    for (auto it = cuts.begin(); it != cuts.end(); it++) {
//        offset_t start, end;
//        if (it->orientation == FRW) {
////            start = left ? it->range.start : 1;
////            end = left ? lr_sketch->sequences[it->seq_id].second - 1 : it->range.end;
//            start = left ? it->range.start : max(0, it->range.start - 700);
//            end = left ? min(lr_sketch->sequences[it->seq_id].second - 1, it->range.end + 700) : it->range.end;
//        } else {
//            start = left ? max(0, it->range.start - 700) : it->range.start;
//            end = left ? it->range.end : min(lr_sketch->sequences[it->seq_id].second - 1, it->range.end + 700);
////            start = left ? 1 : it->range.start;
////            end = left ? it->range.end : lr_sketch->sequences[it->seq_id].second - 1;
//        }
//        it->range.start = start;
//        it->range.end = end;
//        if (it->orientation == FRW)
//            Logger::instance().debug("%s: %d-%d FRW\n", lr_sketch->sequences[it->seq_id].first.c_str(), start, end);
//        else
//            Logger::instance().debug("%s: %d-%d REV\n", lr_sketch->sequences[it->seq_id].first.c_str(), start, end);
//            ans.insert({it->seq_id, *it});
//    }
//
//    return ans;
//}

string InsertionAssembler::get_overlap(vector<string>& l, vector<string>& r, string& m) {
    spoa::Graph graph{};
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kOV, 2, -32, -64, -1);
    graph.Clear();
    for (int i = 0; i < l.size(); i++) {
        auto alignment = alignment_engine->Align(l[i], graph);
        graph.AddAlignment(alignment, l[i]);
    }
    if (m != "") {
        auto alignment = alignment_engine->Align(m, graph);
        graph.AddAlignment(alignment, m);
    }
    for (int i = r.size() - 1; i >= 0; i--) {
        auto alignment = alignment_engine->Align(r[i], graph);
        graph.AddAlignment(alignment, r[i]);
    }

    auto msa = graph.GenerateMultipleSequenceAlignment(true);

    #ifdef DEBUG
    Logger::instance().debug("- MSA\n");
    for (int i = 0; i < msa.size(); i++) {
        Logger::instance().debug("%s\n", msa[i].c_str());
    }
    Logger::instance().debug("\n");
    #endif

    string consensus = msa[msa.size() - 1];
    consensus.erase(std::remove(consensus.begin(), consensus.end(), '-'), consensus.end());

    return consensus;
}

//pair<string, int> InsertionAssembler::assemble(vector<string>& left_reads, vector<string>& right_reads, string &output) {
//	output += "--- Building Long Insertion ---";
//
//    string left_seg, right_seg, lanchor, ranchor, lcheck, rcheck;
//    vector<string> lsegs, rsegs;
//    map<id_t, cut> lcuts, rcuts, mid;
//    unordered_set<hash_t> l_minimizers_frw, l_minimizers_rev, r_minimizers_frw, r_minimizers_rev;
//    int cut_size, support = left_reads.size() + right_reads.size();
//    int step = 1;
//
//    string l_tmp, r_tmp;
//
//    while (true) {
//#ifdef DEBUG
//        Logger::instance().debug("----------------\n");
//        Logger::instance().debug("STEP %d\n", step);
//#endif
//        support += left_reads.size() + right_reads.size();
//
//        left_seg = build_segment(left_reads);
//        lsegs.push_back(left_seg);
//        cut_size = min((int) left_seg.size() / 2, 300);
//        if (step != 0)
//            lanchor = left_seg.substr(left_seg.size() - cut_size, cut_size);
//        else
//            lanchor = left_seg;
//        lcheck = left_seg.substr(0, left_seg.size() - cut_size);
//		lcuts = find_cuts(lanchor, lcheck, true);
//
//        right_seg = build_segment(right_reads);
//        rsegs.push_back(right_seg);
//        cut_size = min((int) right_seg.size() / 2, 300);
//        if (step != 0)
//            ranchor = right_seg.substr(0, cut_size);
//        else
//            ranchor = right_seg;
//        rcheck = right_seg.substr(cut_size, right_seg.size() - cut_size);
//        rcuts = find_cuts(ranchor, rcheck, false);
//
//        output += "\n* " + to_string(step) + ":\n";
//        output += "          Left Segment :  " + left_seg + "\n";
//        output += "          Right Segment:  " + right_seg + "\n";
//        output += "          Picked " + to_string(lcuts.size()) + " new left cuts\n";
//        output += "          Picked " + to_string(rcuts.size()) + " new right cuts\n";
//
//        mid = check_end(lcuts, rcuts);
//        if (mid.size() > 2) {
//			output += "\n===> Insertion built in "+ to_string(step) + " steps.";
//			break;
//		}
//		else if (step == 10) {
//			output += "\n===> Could not build insertion properly in " + to_string(step) + " steps.";
//			break;
//		}
//
//        left_reads = extract_reads(lcuts);
//		Logger::instance().debug("-----\n");
//        right_reads = extract_reads(rcuts);
//
//        step += 1;
//    }
//
//    vector<string> middle_reads = extract_reads(mid);
//    string middle = build_segment(middle_reads);
//    support += middle_reads.size();
//
//#ifdef DEBUG
//    Logger::instance().debug("SEGMENTS:\n");
//    for (int i = 0; i < lsegs.size(); i++)
//        Logger::instance().debug(">l%d\n%s\n", i+1, lsegs[i].c_str());
//    for (int i = 0; i < rsegs.size(); i++)
//        Logger::instance().debug(">r%d\n%s\n", i+1, rsegs[i].c_str());
//    Logger::instance().debug(">m\n%s\n", middle.c_str());
//#endif
//
//    string insertion = get_overlap(lsegs, rsegs, middle);
//
//    return {insertion, support};
//}

pair<string, int> InsertionAssembler::assemble(vector<string>& left_reads, vector<string>& right_reads, string &output) {
    output += "--- Building Long Insertion ---";

    string left_seg, right_seg, lanchor, ranchor, lcheck, rcheck;
    vector<string> lsegs, rsegs;
    unordered_set<hash_t> l_minimizers_frw, l_minimizers_rev, r_minimizers_frw, r_minimizers_rev;
    int cut_size, support = left_reads.size() + right_reads.size(), anchor_cut;
    int step = 1;

    string l_tmp, r_tmp;

    vector<cut> read_cuts;
    vector<cut> lcuts, rcuts, bimodal_cuts;

    bool overlapping_f = false;

    while (true) {
        Logger::instance().debug("STEP %d\n", step);
        if (step == 10) {
            output += "\n===> Could not build insertion properly in " + to_string(step) + " steps.";
            break;
        }

#ifdef DEBUG
        Logger::instance().debug("----------------\n");
        Logger::instance().debug("STEP %d\n", step);
#endif
        support += left_reads.size() + right_reads.size();

        //cerr << "Building Left" << endl;
        left_seg = build_segment(left_reads);
        if (left_seg == "")
            break;
        //cerr << "Finishing Left" << endl;
        lsegs.push_back(left_seg);
        anchor_cut = min(600, (int)left_seg.size());
        //cerr << "Anchor size: " << anchor_cut << endl;
        lanchor = left_seg.substr(left_seg.size() - anchor_cut, anchor_cut - 1);
        cerr << ">s" << step << "-l\n" << lanchor << endl;
        //cerr << "Got Anchor" << endl;

        //cerr << "Building Right" << endl;
        right_seg = build_segment(right_reads);
        if (right_seg == "")
            break;
        //cerr << "Finishing Right" << endl;
        rsegs.push_back(right_seg);
        anchor_cut = min(600, (int)right_seg.size());
        //cerr << "Anchor size: " << anchor_cut << endl;
        ranchor = right_seg.substr(0, anchor_cut - 1);
        cerr << ">s" << step << "-r\n" << ranchor << endl;
        ////cerr << "Got Anchor" << endl;

        unordered_set<hash_t> dummy;
        vector<string> reads;
        pair<vector<cut>, int> cut_results;
        cut_results = lr_sketch->find_cuts_all(lanchor, ranchor, true, reads);
        read_cuts = cut_results.first;

        for (int i = 0; i < read_cuts.size(); i++) {
            if (read_cuts[i].type == PARTIAL_LEFT && lcuts.size() < 7)
                lcuts.push_back(read_cuts[i]);
            else if (read_cuts[i].type == PARTIAL_RIGHT && rcuts.size() < 7)
                rcuts.push_back(read_cuts[i]);
            else if (read_cuts[i].type == BIMODAL && bimodal_cuts.size() < 7)
                bimodal_cuts.push_back(read_cuts[i]);
        }

        if (bimodal_cuts.size() > 2) {
            output += "\n===> Insertion built in "+ to_string(step) + " steps.";
            break;
        }
        else if (cut_results.second == -1) {
            Logger::instance().debug("Overlap\n");
            overlapping_f = true;
            break;
        }

        left_reads = extract_reads(lcuts);
        Logger::instance().debug("-----\n");
        right_reads = extract_reads(rcuts);

        lcuts.clear();
        rcuts.clear();
        bimodal_cuts.clear();

        output += "\n* " + to_string(step) + ":\n";
        output += "          Left Segment :  " + left_seg + "\n";
        output += "          Right Segment:  " + right_seg + "\n";
        output += "          Picked " + to_string(left_reads.size()) + " new left cuts\n";
        output += "          Picked " + to_string(right_reads.size()) + " new right cuts\n";

        step += 1;
    }

    vector<string> middle_reads;
    string middle = "";
    if (!overlapping_f) {
        middle_reads = extract_reads(bimodal_cuts);
        middle = build_segment(middle_reads);
    }
    support += middle_reads.size();

#ifdef DEBUG
    Logger::instance().debug("SEGMENTS:\n");
    for (int i = 0; i < lsegs.size(); i++)
        Logger::instance().debug(">l%d\n%s\n", i+1, lsegs[i].c_str());
    for (int i = 0; i < rsegs.size(); i++)
        Logger::instance().debug(">r%d\n%s\n", i+1, rsegs[i].c_str());
    Logger::instance().debug(">m\n%s\n", middle.c_str());
#endif

    string insertion = get_overlap(lsegs, rsegs, middle);

    return {insertion, support};
}