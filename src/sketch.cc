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
const int MIN_FREQ = 16;
const int MIN_HITS = 48;

Sketch::Sketch(string path, int k, int w) {
	kmer_size = k;
	window_size = w;
	file = path;
    ref_minimizers_vec.reserve(50000);
	build_sketch();
}

void Sketch::sketch_query(std::vector<std::string> reads, int k, int w) {
	kmer_size = k;
	window_size = w;
	build_sketch(reads);
}

void Sketch::build_sketch() {
    int id = 0;
    ifstream fin;
    fin.open(file);
    string tmp;
    int cnt = 0;
    while (!fin.eof()) {
        getline(fin, tmp);
        if (tmp[0] == '>') {
            sequences_names.insert({id, tmp.substr(1, tmp.size() - 1)});
            cnt++;
        } else {
            transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            get_ref_minimizers(&tmp[0], id, tmp.size());
            id++;
        }
    }
    update_ref_sketch();
}

void Sketch::build_sketch(vector<string> reads) {
    for (int i = 0; i < reads.size(); i++) {
        transform(reads[i].begin(), reads[i].end(), reads[i].begin(), ::toupper);
        get_query_minimizers(&reads[i][0], i, reads[i].size());
    }
    update_query_sketch();
    query_minimizers_vec.clear();
}

void Sketch::get_ref_minimizers(char* read, int id, int len) {
    deque<pair<uint64_t, Location> > window;

    char tmp[16];

    int win = window_size - 1;
    for (int i = 0; i < len - kmer_size + 1; i++) {

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

        if (i >= win) {
            if (ref_minimizers_vec.empty() || ref_minimizers_vec.back().first != window.front().first)
                ref_minimizers_vec.push_back(window.front());
        }
    }
}

void Sketch::get_query_minimizers(char* read, int id, int len) {
    deque<pair<uint64_t, int> > window;

    char tmp[16];

    int win = window_size - 1;
    for (int i = 0; i < len - kmer_size + 1; i++) {

        MurmurHash3_x64_128(read + i, kmer_size, 50, tmp);
        uint64_t hash_fw = *((uint64_t *)tmp);

        while(!window.empty() && window.back().first >= hash_fw)
            window.pop_back();

        window.push_back(pair<uint64_t, int>(hash_fw, i));

        while (!window.empty() && window.front().second <= i - window_size)
            window.pop_front();

        if (i >= win) {
            if (query_minimizers_vec.empty() || query_minimizers_vec.back().first != window.front().first)
                query_minimizers_vec.push_back(window.front());
        }
    }
}

void Sketch::update_ref_sketch() {
    sort(ref_minimizers_vec.begin(), ref_minimizers_vec.end(), [](const auto &a, const auto &b) { return a.first < b.first;});
    auto idx = ref_minimizers.begin();
    uint64_t last_hash = UINT64_MAX;
    for (auto it = ref_minimizers_vec.begin(); it != ref_minimizers_vec.end(); it++) {
        if (it->first == last_hash) {
            idx->second.push_back(it->second);
        }
        else {
            last_hash = it->first;
            idx = ref_minimizers.find(it->first);
            if (idx == ref_minimizers.end()) {
                vector<Location> tmp;
                tmp.reserve(128);
                tmp.push_back(it->second);
                idx = ref_minimizers.insert({it->first, tmp}).first;
            }
            else {
                idx->second.push_back(it->second);
            }
        }
    }
}

void Sketch::update_query_sketch() {
    for (int i = 0; i < query_minimizers_vec.size(); i++) {
        query_minimizers.insert(query_minimizers_vec[i].first);
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
    buckets.reserve(BUCKET_SIZE);
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
    freq.reserve(BUCKET_SIZE);
    for (int i = 0; i < (int)buckets.size(); i++) {
        freq.push_back(cnt[i]);
    }

    //first range
    int mx_idx = max_element(freq.begin(), freq.end()) - freq.begin();
    if (freq[mx_idx] < MIN_FREQ) {
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
    pair<int, int> range_1 =  make_pair(l, r);

    //second range
    pair<int, int> range_2;
    mx_idx = max_element(freq.begin() + r, freq.end()) - freq.begin();
    if (freq[mx_idx] < MIN_FREQ) {
        range_2 = make_pair(0, 0);
    }
    else {
        l = mx_idx;
        while (l > 0) {
            if (freq[l - 1] > avg) {
                l--;
            } else {
                break;
            }
        }
        r = mx_idx;
        while (r < (int)freq.size() && freq[r] > avg) {
            r++;
        }
        range_2 =  make_pair(l, r);
    }

    pair<int, int> ans;
    if (buckets[range_1.first].first == 0)
        ans.first = 1;
    else
        ans.first = buckets[range_1.first].first;

    if (range_2.first != 0)
        ans.second = buckets[range_2.second - 1].second;
    else
        ans.second = buckets[range_1.second - 1].second;

    return ans;
}

vector<pair<string, pair<int, int> > > Sketch::find_cuts() {
    map<int, vector<int> > hits_id;
    for (auto it = query_minimizers.begin(); it != query_minimizers.end(); it++) {
        auto idx = ref_minimizers.find(*it);
        if (idx == ref_minimizers.end())
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
        if (it->second.size() < MIN_HITS) {
            continue;
        }
        pair<int, int> a = find_range(it->second);
        if (a.first == 0 && a.second == 0)
            continue;
        cuts.push_back(make_pair(sequences_names[it->first], a));
    }

    return cuts;
}