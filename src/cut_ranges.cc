#include <cstring>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include "cut_ranges.h"

using namespace std;

//TODO merge finding cuts

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

void cut_ranges::extract_reads(string lr_path) {
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
	int length = end - start;
	auto i = reads.find(lr_name);
	if (length > i->second.second.size())
		length = i->second.second.size();
	string cut = i->second.second.substr(start - i->second.first, end - start);
	return cut;
}

//TODO ERROR HANDLING, OPTIMIZE
std::pair<std::string, std::pair<int, int>> cut_consensus(vector<string> alignments) {
    int th = (alignments.size() - 1)/2;
    string consensus = alignments[alignments.size() - 1];
    int left = 0, right = 0;
    int gap_cnt = 0;
    for (int idx = 100; idx > 0; idx = idx/2) {
        for (int i = 0; i < alignments.size() - 2; i++) {
            if (alignments[i][idx] == '-')
                gap_cnt++;
        }
        if (gap_cnt > th) {
            left = idx * 2;
            break;
        }
    }
    gap_cnt = 0;
    for (int idx = 100; idx > 0; idx = idx/2) {
        for (int i = 0; i < alignments.size() - 2; i++) {
            if (alignments[i][alignments[i].size() - idx - 1] == '-')
                gap_cnt++;
        }
        if (gap_cnt > th) {
            right = idx * 2;
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

//TODO multiple clusters
vector<vector<string> > cluster_reads(vector<string> msa, int l, int r) {
    int ncol = msa[0].size() - l - r, nreads = msa.size();
    int profile[ncol];
    memset(profile, 0, sizeof(profile[0]) * ncol);

    for (int i = 0; i < ncol; i++) {
        for (int j = 0; j < nreads; j++) {
            if (msa[j][i + l] == '-')
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

    vector<pair<int, int> > candidates;
    int start = 0, end = 0;
    int f = 0, cnt = 0;
    int up_th = 0.6 * (msa.size() - 1), lw_th = 0.4 * (msa.size() - 1);
    for (int i = 0; i < ncol; i++) {
        if (filtered[i] >= lw_th && filtered[i] <= up_th) {
            if (f == 1) {
                cnt++;
                continue;
            }
            else {
                start = i + l;
                f = 1;
            }
        }
        else {
            if (f == 1) {
                if (cnt > 5) {
                    end = i + l;
                    candidates.push_back({start, end});
                    f = 0;
                    cnt = 0;
                }
                else
                    f = 0;
            }
        }
    }

    pair<int, int> ins;
    int maxx = 0;
    for (int i = 0; i < candidates.size(); i++) {
        int tmp = candidates[i].second - candidates[i].first;
        if (tmp > maxx) {
            ins = candidates[i];
            maxx = tmp;
        }
    }

    vector<vector<string>> clusters;
    vector<string> one, two;

    if (ins.first == 0 && ins.second == 0) {
        one.insert(one.end(), msa.begin(), msa.end() - 1);
        clusters.push_back(one);
        clusters.push_back(two);
        return clusters;
    }

    int th = 0.8 * (ins.second - ins.first);

    for (int i = 0; i < msa.size() - 1; i++) {
        if (count(msa[i].begin() + ins.first, msa[i].begin() + ins.second, '-') >= th) {
            msa[i].erase(std::remove(msa[i].begin(), msa[i].end(), '-'), msa[i].end());
            two.push_back(msa[i]);
        }
        else {
            msa[i].erase(std::remove(msa[i].begin(), msa[i].end(), '-'), msa[i].end());
            one.push_back(msa[i]);
        }
    }

    clusters.push_back(one);
    clusters.push_back(two);

    return clusters;
}