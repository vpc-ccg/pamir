#include "cut_ranges.h"

#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "logger.h"

using namespace std;

//TODO merge finding cuts

cut_ranges::cut_ranges() {}

cut_ranges::cut_ranges(const string &lrPath, const string &dat_path, const string &ranges_path, bool build_index) : lr_path(lrPath) {
    string seq_file = dat_path + "/sequences.dat";
    Logger::instance().info("Loading sequence names from: %s\n", seq_file.c_str());
    ifstream seqin(seq_file, ios::in | ios::binary);
    if (!seqin) {
        Logger::instance().info("Could not load file. Quitting.\n");
        exit(1);
    }

    string name;
    uint16_t len;
    uint8_t name_size;
    while(seqin.read(reinterpret_cast<char*>(&name_size), sizeof(uint8_t))) {
        name.resize(name_size);
        seqin.read((char*)(name.c_str()), name_size);
        seqin.read(reinterpret_cast<char*>(&len), sizeof(uint16_t));
        names.push_back(name);
    }
    seqin.close();

    Logger::instance().info("Loading ranges from: %s\n", ranges_path.c_str());
    ifstream ranges_file(ranges_path, ios::in | ios::binary);
    if (!ranges_file) {
        Logger::instance().info("Could not load file. Quitting.\n");
        exit(1);
    }

    uint32_t size;
    ranges_file.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    ranges.resize(size);
    ranges_file.read(reinterpret_cast<char*>(&ranges[0]), size * sizeof(pair<id_t, range_s>));
    ranges_file.close();
}

cut_ranges::cut_ranges(const string &lrPath, bool build_index) : lr_path(lrPath) {
    if (build_index) {
        string lr_name;
        lr_file.open(lr_path);
        string line;

        while (!lr_file.eof()) {
            getline(lr_file, line);
            if (line[0] == '>') {
                lr_name = line.substr(1);
                read_offsets.insert({lr_name, lr_file.tellg()});
//                //cerr << lr_name <<  " | " << lr_file.tellg() << endl;
                getline(lr_file, line);
            }
        }
    }
//    //cerr << "----------------------" << endl;
    lr_file.close();
    lr_file.open(lr_path);
}

void cut_ranges::extract() {
    ifstream fin;
    fin.open(lr_path);
    string line;

    int curr_read = 0;
    vector<pair<id_t, range_s> >::iterator search_read = ranges.begin();

    while (!fin.eof() && search_read != ranges.end()) {
        getline(fin, line);
        getline(fin, line);
        if (curr_read == (*search_read).first) {
            string read = line.substr((*search_read).second.start, (*search_read).second.end);
            reads.push_back({curr_read, {(*search_read).second.start, read}});
            search_read++;
        }
        curr_read++;
    }
}

read_cut cut_ranges::find_read(id_t id) {
    read_cut dummy;
    dummy.first = id;
    vector<read_cut>::iterator it = std::lower_bound(reads.begin(), reads.end(), dummy);
    return (*it);
}

pair<string, string> cut_ranges::get_cut(id_t lr_id, offset_t start, offset_t end) {
    string cut;
    read_cut read = find_read(lr_id);
    cut = read.second.second.substr(start - read.second.first, end - start);
	return {names[lr_id], cut};
}

string cut_ranges::get_cut(string name, offset_t start, offset_t end) {
    string cut;
    lr_file.seekg(read_offsets[name]);
    string line;
    getline(lr_file, line);
    Logger::instance().debug(">%s\n%s\n", name.c_str(), line.c_str());
//    //cerr << read_offsets[name] << endl;
//    //cerr << line << endl;
    cut = line.substr(start, end - start - 1);
    return cut;
}

//std::pair<std::string, std::pair<int, int>> binary_search(int begin, int end, int search_val) {
//    int mid = (begin + end)/2;
//    if (mid < search_val) {
//        binary_search(mid + 1, end, search_val);
//    }
//    else {
//        if (mid - 1 < search_val)
//            return;
//        binary_search(begin, mid - 1);
//    }
//}

//TODO: Check for possible unhandled cases
std::pair<std::string, std::pair<int, int>> cut_consensus_bimodal(vector<string> alignments, int left_reads,
                                                                  int right_reads, int bimodal_reads) {
    int th_left = alignments.size() - 1 - 0.5 * (bimodal_reads + left_reads);
    int th_right = alignments.size() - 1 - 0.5 * (bimodal_reads + right_reads);
    string consensus = alignments[alignments.size() - 1];

    bool up = false, down = false;
    int step = 200;
    int prv = 0;
    int left = 0, right = 0;
    int gap_cnt;
    int curr_answer = 0;
    cerr << "left threshold: " << th_left << " | right threshold: " << th_right << endl << endl;
    cerr << " ----- Left ----- " << endl;
    for (int idx = min(200, (int)alignments[0].size()/2); idx < alignments[0].size(); idx += step) {
        gap_cnt = 0;
        for (int i = 0; i < alignments.size() - 1; i++) {
            if (alignments[i][idx] == '-')
                gap_cnt++;
        }
        cerr << "......" << endl;
        cerr << "idx: " << idx << endl;
        cerr << "gap cnt: " << gap_cnt << endl;
        if (gap_cnt > th_left) {
            if (down) {
                cerr << "curr answer: " << curr_answer << endl;
                if (step != 1)
                    step = abs(curr_answer - idx)/2;
                if (step == 0)
                    step = 1;
                cerr << "down==true" << endl;
                cerr << "new step: " << step << endl;
            }
            else {
                cerr << "curr answer: " << curr_answer << endl;
                if (step != 1)
                    step = abs(curr_answer - idx)/2;
                if (step == 0)
                    step = 1;
//                step = idx;
                cerr << "down==false" << endl;
                cerr << "new step: " << step << endl;
            }
            prv = idx;
            up = true;
            down = false;
            cerr << "up true, down false" << endl;
        }
        else {
            down = true;
            cerr << "down true" << endl;
            if (step != 0 && step != 1) {
                cerr << "up==false" << endl;
                step = -abs(prv - idx) / 2;
                cerr << "new step: " << step << endl;
                prv = idx;
                curr_answer = idx;
            } else {
                cerr << "up==true" << endl;
                cerr << "left set to " << left << endl;
                left = idx;
                break;
            }
        }
    }

    down = false;
    up = false;
    step = 200;
    prv = 0;
    curr_answer = 0;
    cerr << endl << " ----- Right ----- " << endl;
    for (int idx = min(200, (int)alignments[0].size()/2); idx < alignments[0].size(); idx += step) {
        gap_cnt = 0;
        for (int i = 0; i < alignments.size() - 1; i++) {
            if (alignments[i][alignments[i].size() - idx - 1] == '-')
                gap_cnt++;
        }
        cerr << "......" << endl;
        cerr << "idx: " << idx << endl;
        cerr << "gap cnt: " << gap_cnt << endl;
        if (gap_cnt > th_right) {
            if (down) {
                cerr << "curr answer: " << curr_answer << endl;
                if (step != 1)
                    step = abs(curr_answer - idx)/2;
                if (step == 0)
                    step = 1;
                cerr << "down==true" << endl;
                cerr << "new step: " << step << endl;
            }
            else {
                cerr << "curr answer: " << curr_answer << endl;
                if (step != 1)
                    step = abs(curr_answer - idx)/2;
                if (step == 0)
                    step = 1;
//                step = idx;
                cerr << "down==false" << endl;
                cerr << "new step: " << step << endl;
            }
            prv = idx;
            up = true;
            down = false;
            cerr << "up false, down true" << endl;
        }
        else {
            down = true;
            cerr << "down true" << endl;
            if (step != 0 && step != 1) {
                cerr << "up==false" << endl;
                step = -abs(prv - idx) / 2;
                cerr << "new step: " << step << endl;
                prv = idx;
                curr_answer = idx;
            } else {
                cerr << "up==true" << endl;
                right = idx;
                cerr << "right set to " << left << endl;
                break;
            }
        }
    }

    cerr << "Final answer: " << left << ", " << right << endl;
    pair<int, int> ans = {left, right};
    string cut = consensus.substr(left, consensus.size() - right - left);
    cut.erase(std::remove(cut.begin(), cut.end(), '-'), cut.end());
    return {cut, ans};
}

std::pair<std::string, std::pair<int, int>> cut_consensus_single(vector<string> alignments) {
    int th = alignments.size() - 1 - 0.5 * (alignments.size() - 1);
    string consensus = alignments[alignments.size() - 1];

    bool up = false, down = false;
    int left = 0, right = 0;
    int gap_cnt;
    int step = 200, prv = 0;
    int curr_answer = 0;
    //cerr << "threshold: " << th << endl;
    //cerr << " ----- Left ----- " << endl;
    for (int idx = min(200, (int)alignments[0].size()/2); idx < alignments[0].size(); idx += step) {
        gap_cnt = 0;
        for (int i = 0; i < alignments.size() - 1; i++) {
            if (alignments[i][idx] == '-')
                gap_cnt++;
        }
        cerr << "......" << endl;
        cerr << "idx: " << idx << endl;
        cerr << "gap cnt: " << gap_cnt << endl;
        if (gap_cnt > th) {
            if (down) {
                cerr << "curr answer: " << curr_answer << endl;
                if (step != 1)
                    step = abs(curr_answer - idx)/2;
                if (step == 0)
                    step = 1;
                cerr << "down==true" << endl;
                cerr << "new step: " << step << endl;
            }
            else {
                cerr << "curr answer: " << curr_answer << endl;
                if (step != 1)
                    step = abs(curr_answer - idx)/2;
                if (step == 0)
                    step = 1;
//                step = idx;
                cerr << "down==false" << endl;
                cerr << "new step: " << step << endl;
            }
            prv = idx;
            up = true;
            down = false;
            cerr << "up true, down false" << endl;
        }
        else {
            down = true;
            cerr << "down true" << endl;
            if (step != 0 && step != 1) {
                cerr << "up==false" << endl;
                step = -abs(prv - idx) / 2;
                cerr << "new step: " << step << endl;
                prv = idx;
                curr_answer = idx;
            } else {
                cerr << "up==true" << endl;
                cerr << "left set to " << left << endl;
                left = idx;
                break;
            }
        }
    }

    step = 200;
    prv = 0;
    up = false, down = false;
    curr_answer = 0;
    for (int idx = min(200, (int)alignments[0].size()/2); idx < alignments[0].size(); idx += step) {
        gap_cnt = 0;
        for (int i = 0; i < alignments.size() - 1; i++) {
            if (alignments[i][alignments[i].size() - idx - 1] == '-')
                gap_cnt++;
        }
        cerr << "......" << endl;
        cerr << "idx: " << idx << endl;
        cerr << "gap cnt: " << gap_cnt << endl;
        if (gap_cnt > th) {
            if (down) {
                cerr << "curr answer: " << curr_answer << endl;
                if (step != 1)
                    step = abs(curr_answer - idx)/2;
                if (step == 0)
                    step = 1;
                cerr << "down==true" << endl;
                cerr << "new step: " << step << endl;
            }
            else {
                cerr << "curr answer: " << curr_answer << endl;
                if (step != 1)
                    step = abs(curr_answer - idx)/2;
                if (step == 0)
                    step = 1;
//                step = idx;
                cerr << "down==false" << endl;
                cerr << "new step: " << step << endl;
            }
            prv = idx;
            up = true;
            down = false;
            cerr << "up true, down false" << endl;
        }
        else {
            down = true;
            cerr << "down true" << endl;
            if (step != 0 && step != 1) {
                cerr << "up==false" << endl;
                step = -abs(prv - idx) / 2;
                cerr << "new step: " << step << endl;
                prv = idx;
                curr_answer = idx;
            } else {
                cerr << "up==true" << endl;
                cerr << "left set to " << left << endl;
                right = idx;
                break;
            }
        }
    }

    pair<int, int> ans = {left, right};
    string cut = consensus.substr(left, consensus.size() - right - left);
    cut.erase(std::remove(cut.begin(), cut.end(), '-'), cut.end());
    return {cut, ans};
}

void insertionSort(int arr[], int n) {
    int key, j;
    for (int i = 1; i < n; i++)
    {
        key = arr[i];
        j = i - 1;

        while (j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}

//TODO multiple cluster
//TODO merge with cut_consensus
vector<vector<string> > cluster_reads(vector<string> msa, int l, int r) {
    int ncol = msa[0].size();
    int nreads = msa.size() - 1;
    int profile[ncol];
    memset(profile, 0, sizeof(profile[0]) * ncol);

    for (int i = 0; i < ncol; i++) {
        for (int j = 0; j < nreads; j++) {
            if (msa[j][i] == '-')
                profile[i]++;
        }
    }

    int window[5], row = 0, col = 0;
    int filtered[ncol];
    memset(filtered, 0, sizeof(filtered[0]) * ncol);

    for(col = 2; col < ncol - 2; col++) {

        window[0] = profile[col-2];
        window[1] = profile[col-1];
        window[2] = profile[col];
        window[3] = profile[col+1];
        window[4] = profile[col+2];

        insertionSort(window,5);

        filtered[col]=window[2];
    }

    int val_th = 0.3 * nreads;
    int i = 0, val_s, val_e;
    while (profile[i] > val_th)
        i++;
    val_s = i;
    i = ncol - 1;
    while (profile[i] > val_th)
        i--;
    val_e = i;

    vector<pair<int, int> > candidates;
    int start = 0, end = 0;
    int f = 0, cnt = 0;
    int up_th = 0.6 * nreads, lw_th = 0.4 * nreads;
    for (int i = val_s; i <= val_e; i++) {
        if (filtered[i] >= lw_th && filtered[i] <= up_th) {
            if (f == 1) {
                cnt++;
                continue;
            }
            else {
                start = i;
                f = 1;
            }
        }
        else {
            if (f == 1) {
                if (cnt > 5) {
                    end = i;
                    candidates.push_back({start, end});
                    f = 0;
                    cnt = 0;
                }
                else
                    f = 0;
            }
        }
    }

    vector<vector<string>> clusters;
    vector<string> one, two;
    pair<int, int> ins;

    if (candidates.size() == 0) {
        one.insert(one.end(), msa.begin(), msa.end() - 1);
        clusters.push_back(one);
        clusters.push_back(two);
        return clusters;
    }
    else if (candidates.size() == 1) {
        ins = candidates[0];
    }
    else {
        int maxx = 0;
        for (int i = 0; i < candidates.size(); i++) {
            int tmp = candidates[i].second - candidates[i].first;
            if (tmp > maxx) {
                ins = candidates[i];
                maxx = tmp;
            }
        }
    }

    int th = 0.9 * (ins.second - ins.first);

    for (int i = 0; i < nreads; i++) {
        if (count(msa[i].begin(), msa[i].begin() + ins.first, '-') == ins.first) {
            msa[i].erase(std::remove(msa[i].begin(), msa[i].end(), '-'), msa[i].end());
            one.push_back(msa[i]);
        }
        else if (count(msa[i].begin() + ins.first, msa[i].end(), '-') == msa[i].size() - ins.first) {
            msa[i].erase(std::remove(msa[i].begin(), msa[i].end(), '-'), msa[i].end());
            one.push_back(msa[i]);
        }
        else if (count(msa[i].begin() + ins.first, msa[i].begin() + ins.second, '-') >= th) {
            msa[i].erase(std::remove(msa[i].begin(), msa[i].end(), '-'), msa[i].end());
            two.push_back(msa[i]);
        }
        else {
            msa[i].erase(std::remove(msa[i].begin(), msa[i].end(), '-'), msa[i].end());
            one.push_back(msa[i]);
        }
    }

    if (!two.empty() && two.size() == 1) {
        one.push_back(two[0]);
        two.clear();
        clusters.push_back(one);
        clusters.push_back(two);
    }
    else {
        clusters.push_back(one);
        clusters.push_back(two);
    }

    return clusters;
}
