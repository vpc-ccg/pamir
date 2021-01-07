#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "cut_ranges.h"

using namespace std;

//TODO merge finding cuts

cut_ranges::cut_ranges() {}

cut_ranges::cut_ranges(const string &lrPath, bool build_index) : lr_path(lrPath) {
    if (build_index) {
        string lr_name;
        ifstream fin;
        fin.open(lr_path);
        string line;

        while (!fin.eof()) {
            getline(fin, line);
            if (line[0] == '>') {
                lr_name = line.substr(1, line.size() - 1);
                read_offsets.insert({lr_name, fin.tellg()});
                getline(fin, line);
            }
        }
    }
}

void cut_ranges::add_range(string lr, int start, int end) {
	auto curr_lr = lr_ranges.find(lr);
	if (curr_lr == lr_ranges.end())
		lr_ranges.insert(pair<string, pair<int, int>>(lr, pair<int, int>(start, end)));
	else {
		if (start < curr_lr->second.first)
			curr_lr->second.first = start;
		if (end > curr_lr->second.second)
			curr_lr->second.second = end;
	}
}

void cut_ranges::extract() {
	string lr_name;
	ifstream fin;
	fin.open(lr_path);
	string line;

	while (!fin.eof()) {
		getline(fin, line);
		if (line[0] == '>') {
		    //TODO fix
			istringstream iss(line);
    		vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
			lr_name = tokens[0].substr(1, line.size() - 1);
			auto i = lr_ranges.find(lr_name);
			if (i != lr_ranges.end()) {
				getline(fin, line);
				int end = i->second.second;
				if (end > line.size())
					end = line.size();
				string read = line.substr(i->second.first - 1, end - i->second.first);
				reads.insert(pair<string, pair<int, string> >(lr_name, pair<int, string> (i->second.first, read)));
			}
		}
	}
}

string cut_ranges::get_cut(string lr_name, int start, int end) {
    string cut;
    if (read_offsets.empty()) {
        int length = end - start;
        auto i = reads.find(lr_name);
        if (length > i->second.second.size())
            length = i->second.second.size();
        cut = i->second.second.substr(start - i->second.first, end - start);
    }
	else {
        ifstream fin;
        fin.open(lr_path);
        fin.seekg(read_offsets[lr_name]);
        string line;
        getline(fin, line);
		int s = start;
		if (s == 1)
			s = 0;
        cut = line.substr(s, end - s);
	}
	return cut;
}

//TODO ERROR HANDLING, OPTIMIZE
std::pair<std::string, std::pair<int, int>> cut_consensus_bimodal(vector<string> alignments, int left_reads,
		int right_reads, int bimodal_reads) {
    int cnt = alignments.size() - 1;
    int th_left = 0.5 * (bimodal_reads + left_reads);
    int th_right = 0.5 * (bimodal_reads + right_reads);
    string consensus = alignments[alignments.size() - 1];

    int left = 0, right = 0;
    int gap_cnt;
    for (int idx = 0; idx < alignments[0].size(); idx = idx + 50) {
        gap_cnt = 0;
        for (int i = 0; i < alignments.size() - 1; i++) {
            if (alignments[i][idx] == '-')
                gap_cnt++;
        }
        if (gap_cnt < th_left) {
            left = idx;
            break;
        }
    }

    for (int idx = 0; idx < alignments[0].size(); idx = idx + 50) {
        gap_cnt = 0;
        for (int i = 0; i < alignments.size() - 1; i++) {
            if (alignments[i][alignments[i].size() - idx - 1] == '-')
                gap_cnt++;
        }
        if (gap_cnt < th_right) {
            right = idx;
            break;
        }
    }
    pair<int, int> ans = {left, right};
    string cut = consensus.substr(left, consensus.size() - right - left);
    cut.erase(std::remove(cut.begin(), cut.end(), '-'), cut.end());
    return {cut, ans};
}

std::pair<std::string, std::pair<int, int>> cut_consensus_single(vector<string> alignments) {
    int th = 0.6 * (alignments.size() - 1);
    string consensus = alignments[alignments.size() - 1];
    int left = 0, right = 0;
    int gap_cnt;
    for (int idx = 0; idx < alignments[0].size(); idx = idx + 50) {
        gap_cnt = 0;
        for (int i = 0; i < alignments.size() - 1; i++) {
            if (alignments[i][idx] == '-')
                gap_cnt++;
        }
        if (gap_cnt < th) {
            left = idx;
            break;
        }
    }

    for (int idx = 0; idx < alignments[0].size(); idx = idx + 50) {
        gap_cnt = 0;
        for (int i = 0; i < alignments.size() - 1; i++) {
            if (alignments[i][alignments[i].size() - idx - 1] == '-')
                gap_cnt++;
        }
        if (gap_cnt < th) {
            right = idx;
            break;
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
