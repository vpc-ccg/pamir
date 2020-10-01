#include <deque>
#include <string>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "MurmurHash3.h"
#include "sketch.h"

using namespace std;

//TODO FIX
const int MAXN = 30000;
const int BUCKET_SIZE = 100;
const int TH = 90;

Sketch::Sketch(string path, int k, int w) {
	kmer_size = k;
	window_size = w;
	file = path;
	build_sketch();
}

Sketch::Sketch(vector<string> reads, int k, int w) : file(NULL) {
	kmer_size = k;
	window_size = w;
	build_sketch(reads);
}

void Sketch::build_sketch() {
	int id = 0;
	ifstream fin;
	fin.open(file);
	string tmp;
	while (!fin.eof()) {
		getline(fin, tmp);
		if (tmp[0] == '>') {
			sequences.insert(pair<int, string>(id, tmp.substr(1, tmp.size() - 1)));
		}
		else {
			transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
			get_minimizers(&tmp[0], id, tmp.size());
			id++;
		}
	}
}

void Sketch::build_sketch(vector<string> reads) {
	for (int i = 0; i < reads.size(); i++) {
		transform(reads[i].begin(), reads[i].end(), reads[i].begin(), ::toupper);
		get_minimizers(&reads[i][0], -1, reads[i].size());
	}
}

void Sketch::get_minimizers(char* read, int id, int len) {
	deque<pair<uint64_t, Location> > window;

	for (int i = 0; i < len - kmer_size + 1; i++) {
		char tmp[16];
		MurmurHash3_x64_128(read + i, kmer_size, 50, tmp);
		uint64_t hash_fw = *((uint64_t *)tmp);

		Location loc;
		loc.seq_id = id;
		loc.offset = i;

		while(!window.empty() && window.back().first >= hash_fw)
			window.pop_back();

		window.push_back(pair<uint64_t, Location>(hash_fw, loc));

		while (!window.empty() && window.front().second.offset <= i - window_size)
			window.pop_front();

		if (i >= window_size - 1) {
			pair<uint64_t, Location> mini = window.front();
			auto j = minimizers.find(mini.first);
			if (j == minimizers.end()) {
				vector<Location> tmp;
				tmp.push_back(mini.second);
				minimizers.insert(pair<uint64_t, vector<Location> >(mini.first, tmp));				
			}
			else {
				if (!(j->second.back().seq_id == mini.second.seq_id && j->second.back().offset == i))
					j->second.push_back(mini.second);
			}
		}
	}
}

void Sketch::print() {
	for(auto it = minimizers.begin(); it != minimizers.end(); it++) {
    	cout << it->first << "\t";
		for (int i = 0; i < it->second.size(); i++) {
			cout << "\t" << it->second[i].offset << "\t" << it->second[i].seq_id << "\n"; 
		}
	}
}

int get_avg(vector<int> v) {
	int sum = 0;
	for (int i = 0; i < (int)v.size(); i++) {
		sum += v[i];
	}
	int res = sum / v.size();
	return res;
}

pair<int, int> find_range(vector<int> hits) {
	sort(hits.begin(), hits.end());
	vector<pair<int, int> > buckets;
	map<int, int> cnt;
	for (int s = 0; s < MAXN; s += BUCKET_SIZE) {
		buckets.push_back(make_pair(s, s + BUCKET_SIZE));
	}
	for (int i = 0, j = 0; i < (int)hits.size(); i++) {
		while (hits[i] >= buckets[j].second) {
			j++;
		}
		assert(hits[i] >= buckets[j].first);
		cnt[j]++;
	}	
	vector<int> freq;
	for (int i = 0; i < (int)buckets.size(); i++) {
		freq.push_back(cnt[i]);
	}
	int mx_idx = max_element(freq.begin(), freq.end()) - freq.begin();
	
	if (freq[mx_idx] < TH) {
		return make_pair(0, 0);
	}

	int avg = get_avg(freq);

	int l = mx_idx;
	while (l > 0) {
		if (freq[l - 1] > avg) {
			l--;
		} else {
			break;
		}
	}
	int r = mx_idx;
	while (r < (int)freq.size() && freq[r] > avg) {
		r++;
	}
	pair<int, int> range =  make_pair(l, r);
	pair<int, int> ans;
	if (buckets[range.first].first == 0) {
		ans.first = 1;
		ans.second = buckets[range.second - 1].second;
	}
	else {
		ans.first = buckets[range.first].first;
		ans.second = buckets[range.second - 1].second;
	}
	return ans;
}

vector<pair<string, pair<int, int> > > Sketch::find_cuts(Sketch q) {
	map<int, vector<int> > hits_id;
	for (auto it = q.minimizers.begin(); it != q.minimizers.end(); it++) {
		auto idx = minimizers.find(it->first);
		if (idx == minimizers.end())
			continue;
		for (auto i = idx->second.begin(); i != idx->second.end(); i++) {
			auto j = hits_id.find(i->seq_id);
			if (j == hits_id.end()) {
				vector<int> tmp;
				tmp.reserve(32);
				tmp.push_back(i->offset);
				hits_id.insert(pair<int, vector<int> >(i->seq_id, tmp));				
			}
			else {
				j->second.push_back(i->offset);
			}
		}
	}
	
	vector<pair<string, pair<int, int> > > cuts;
	for (auto it = hits_id.begin(); it != hits_id.end(); it++) {
		if (it->second.size() < 300) {
			continue;
		}
		pair<int, int> a = find_range(it->second);
		if (a.first == 0 && a.second == 0)
			continue;
		cuts.push_back(make_pair(sequences[it->first], a));
	}

	return cuts;
}