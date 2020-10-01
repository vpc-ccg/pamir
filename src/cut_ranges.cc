#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include "cut_ranges.h"

using namespace std;

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

string cut_consensus(vector<string> alignments) {
	int th = (alignments.size() - 1)/2 + 1;
	string consensus = alignments[alignments.size() - 1];
	int left = 0, right = consensus.size();
	int gap_cnt = 0;
	for (int idx = 100; idx >= 0; idx = idx/2) {
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
	for (int idx = 100; idx >= 0; idx = idx/2) {
		for (int i = 0; i < alignments.size() - 2; i++) {
			if (alignments[i][alignments[i].size() - idx - 1] == '-')
				gap_cnt++; 
		}
		if (gap_cnt > th) {
			right = idx * 2;
			break;
		}
	}
	string cut = consensus.substr(left, consensus.size() - right - left);
	cut.erase(std::remove(cut.begin(), cut.end(), '-'), cut.end());
	return cut;
}