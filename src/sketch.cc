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

Sketch::Sketch() {}

Sketch::Sketch(string path, bool l, int k, int w) {
    kmer_size = k;
    window_size = w;
    file_path = path;
    if (!l) {
        ref_minimizers_vec.reserve(50000);
        build_sketch();
        dump();
    }
    else {
        load();
    }
}

void Sketch::dump() {
    ofstream fout("sketch.dat", ios::out | ios::binary);
    for (auto it = ref_minimizers.begin(); it != ref_minimizers.end(); it++) {
        typename vector<Location>::size_type size = it->second.size();
        fout.write((char*)&it->first, sizeof(uint64_t));
        fout.write((char*)&size, sizeof(size));
        fout.write((char*)&(it->second[0]), size * sizeof(Location));
    }
    fout.close();
    ofstream seqout("sequences.dat", ios::out | ios::binary);
    for (auto it = sequences.begin(); it != sequences.end(); it++) {
        uint8_t size = it->second.first.length();
        seqout.write((char*)&it->first, sizeof(uint32_t));
        seqout.write((char*)&size, sizeof(uint8_t));
        seqout.write(it->second.first.c_str(), size);
        seqout.write((char*)&it->second.second, sizeof(uint16_t));
    }
    seqout.close();
}

void Sketch::load() {
    ifstream fin(file_path + "/sketch.dat", ios::in | ios::binary);
    uint64_t hash;
    while(fin.read(reinterpret_cast<char*>(&hash), sizeof(uint64_t))) {
        typename vector<Location>::size_type size = 0;
        fin.read(reinterpret_cast<char*>(&size), sizeof(size));
        vector<Location> tmp;
        tmp.resize(size);
        fin.read(reinterpret_cast<char*>(&tmp[0]), size * sizeof(Location));
        ref_minimizers.insert({hash, tmp});
    }
    fin.close();
    ifstream seqin(file_path + "/sequences.dat", ios::in | ios::binary);
    uint32_t id;
    string name;
    uint16_t len;
    while(seqin.read(reinterpret_cast<char*>(&id), sizeof(uint32_t))) {
        uint8_t size;
        seqin.read(reinterpret_cast<char*>(&size), sizeof(uint8_t));
        name.resize(size);
        seqin.read((char*)(name.c_str()), size);
        seqin.read(reinterpret_cast<char*>(&len), sizeof(uint16_t));
        sequences.insert({id, {name, len}});
    }
    seqin.close();
}

void Sketch::sketch_query(vector<string> reads, int k, int w) {
    query_minimizers.clear();
	kmer_size = k;
	window_size = w;
	build_sketch(reads);
}

void Sketch::sketch_query(string query, int k, int w) {
    query_minimizers.clear();
    kmer_size = k;
    window_size = w;
    vector<string> queries;
    queries.push_back(query);
    build_sketch(queries);
}

void Sketch::build_sketch() {
    int id = 0;
    ifstream fin;
    fin.open(file_path);
    string tmp;
    int cnt = 0;
    string name;
    int size;
    while (!fin.eof()) {
        getline(fin, tmp);
        if (tmp[0] == '>') {
            name = tmp.substr(1, tmp.size() - 1);
            cnt++;
        } else {
            size = tmp.size();
            sequences.insert({id, {name, size}});
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

cut find_range(vector<pair<int, uint64_t> > hits) {
    sort(hits.begin(), hits.end(), [](auto &a, auto &b) {
        return a.first < b.first;
    });

    //Create buckets
    vector<pair<int, int> > buckets;
    buckets.reserve(BUCKET_SIZE);
    map<int, int> cnt;
    for (int s = 0; s < MAXN; s += BUCKET_SIZE) {
        buckets.push_back(make_pair(s, s + BUCKET_SIZE));
    }

    //Populate buckets by hit locations
    for (int i = 0, j = 0; i < (int)hits.size(); i++) {
        while (hits[i].first >= buckets[j].second) {
            j++;
        }
        assert(hits[i].first >= buckets[j].first);
        cnt[j]++;
    }

    //Find bucket weights
    vector<int> freq;
    freq.reserve(BUCKET_SIZE);
    for (int i = 0; i < (int)buckets.size(); i++) {
        freq.push_back(cnt[i]);
    }

    //find cut around highest peak
    int mx_idx = max_element(freq.begin(), freq.end()) - freq.begin();

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


    //find cut around second highest peak
    pair<int, int> range_2;
    //find if peak is on the left side or the right side of maximum peak
    int mx_idx_l = max_element(freq.begin(), freq.begin() + max(l - 1, 0)) - freq.begin();
    int mx_idx_r = max_element(freq.begin() + min(r + 1, int (freq.size())), freq.end()) - freq.begin();

    int prv_peak = freq[mx_idx];
    mx_idx = freq[mx_idx_l] > freq[mx_idx_r] ? mx_idx_l : mx_idx_r;

    if (mx_idx <= r && mx_idx_l >= l)
        range_2 = make_pair(0, 0);
    //Discard if new peak is too low comparing to the maximum peak (probably doesn't belong to this insertion)
    else if (freq[mx_idx] < prv_peak/2) {
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

    pair<int, int> final_range;
    int type;
    cut ans;

    //Case 1: only one peak found
    if (range_2.first == range_2.second) {
        if (buckets[range_1.first].first == 0)
            final_range.first = 1;
        else
            final_range.first = buckets[range_1.first].first;
        final_range.second = buckets[range_1.second - 1].second;
        type = SINGLE_PEAK;
        ans.peak1 = final_range;
        ans.peak2 = {0, 0};
    }
    //Case 2: Two peaks found/Bimodal read
    else {
        type = BIMODAL;
        //If second peak is on the right side of the first peak
        if (range_1.first < range_2.first) {
            if (buckets[range_1.first].first == 0)
                final_range.first = 1;
            else
                final_range.first = buckets[range_1.first].first;
            final_range.second = buckets[range_2.second - 1].second;
            ans.peak1 = {buckets[range_1.first].first, buckets[range_1.second - 1].second};
            ans.peak2 = {buckets[range_2.first].first, buckets[range_2.second - 1].second};
        }
        //If second peak is on the left side of the first peak
        else {
            if (buckets[range_2.first].first == 0)
                final_range.first = 1;
            else
                final_range.first = buckets[range_2.first].first;
            final_range.second = buckets[range_1.second - 1].second;
            ans.peak2 = {buckets[range_1.first].first, buckets[range_1.second - 1].second};
            ans.peak1 = {buckets[range_2.first].first, buckets[range_2.second - 1].second};
        }

    }

    ans.range = final_range;
    ans.type = type;

    return ans;
}

float minimizer_similarity(unordered_set<uint64_t> ref, vector<pair<int, uint64_t> > q) {
    int cnt = 0;
    for (int i = 0; i < q.size(); i++) {
        if (find(ref.begin(), ref.end(), q[i].second) != ref.end())
            cnt += 1;
    }
    return float(cnt)/float(q.size());
}

vector<pair<string, pair<pair<int, int>, int> > > Sketch::classify_reads(map<int, vector<pair<int, uint64_t> > > hits,
                                                                         vector<pair<int, cut> > cuts_tmp) {
    //Bimodal reads first, then ordered by number of minimizers used to pick the reads
    sort(cuts_tmp.begin(), cuts_tmp.end(), [](auto &a, auto &b) {
        if (a.second.type < b.second.type)
            return true;
        else if (a.second.type == b.second.type)
            return a.second.number_of_minimizers > b.second.number_of_minimizers;
        else
            return false;
    });

    unordered_set<uint64_t> left_minimizers, right_minimizers;
    vector<pair<string, pair<pair<int, int>, int> > > cuts;

    int i = 0;
    //Separate left and right minimizers
    while (cuts_tmp[i].second.type == BIMODAL) {
        int curr_read = cuts_tmp[i].first;
        cut curr_cut = cuts_tmp[i].second;

        for (auto it = hits[curr_read].begin(); it != hits[curr_read].end(); it++) {
            //Left minimizers
            if (it->first < curr_cut.peak1.second)
                left_minimizers.insert(it->second);
                //Right minimizers
            else if (it->first > curr_cut.peak2.first)
                right_minimizers.insert(it->second);
        }
        cuts.push_back({sequences[curr_read].first, {curr_cut.range, BIMODAL}});
        i += 1;
    }

    for (i; i < cuts_tmp.size(); i++) {
        int curr_read = cuts_tmp[i].first;
        cut curr_cut = cuts_tmp[i].second;
        //Case 1: We have seen a bimodal read or we have already populated the left and right sets
        // -> just compare similarity to left and right sets
        if (!left_minimizers.empty() && !right_minimizers.empty()) {
            float left_similarity = minimizer_similarity(left_minimizers, hits[curr_read]);
            float right_similarity = minimizer_similarity(right_minimizers, hits[curr_read]);

            if (left_similarity > right_similarity)
                cuts.push_back({sequences[curr_read].first, {curr_cut.range, LEFT}});
            else if (right_similarity > left_similarity)
                cuts.push_back({sequences[curr_read].first, {curr_cut.range, RIGHT}});
            else
                cuts.push_back({sequences[curr_read].first, {curr_cut.range, MISC}});
        }
            //Case 2: No bimodal reads, must classify reads on the go -> first set of seen minimizers are assigned to left
        else if (left_minimizers.empty()) {
            for (int j = 0; j < hits[curr_read].size(); j++) {
                left_minimizers.insert(hits[curr_read][j].second);
            }
            cuts.push_back({sequences[curr_read].first, {curr_cut.range, LEFT}});
        }
            //Case 3: all seen reads have been classified as left and right is still empty -> if this new read doesn't have
            //        enough similarity with the left set, assign it to right
        else if (right_minimizers.empty()) {
            float left_similarity = minimizer_similarity(left_minimizers, hits[curr_read]);
            if (left_similarity > 0.5) {
                cuts.push_back({sequences[curr_read].first, {curr_cut.range, LEFT}});
                continue;
            }
            for (int j = 0; j < hits[curr_read].size(); j++) {
                right_minimizers.insert(hits[curr_read][j].second);
            }
            cuts.push_back({sequences[curr_read].first, {curr_cut.range, RIGHT}});
        }
    }

    return cuts;
}

vector<pair<string, pair<pair<int, int>, int> > > Sketch::find_cuts(bool classify) {
    map<int, vector<pair<int, uint64_t> > > hits;

    //Find all hits
    for (auto it = query_minimizers.begin(); it != query_minimizers.end(); it++) {
        auto idx = ref_minimizers.find(*it);
        if (idx == ref_minimizers.end())
            continue;
        for (auto i = idx->second.begin(); i != idx->second.end(); i++) {
            auto j = hits.find(i->seq_id);
            if (j == hits.end()) {
                vector<pair<int, uint64_t> > tmp;
                tmp.reserve(32);
                tmp.push_back({i->offset, *it});
                hits.insert({i->seq_id, tmp});
            }
            else {
                j->second.push_back({i->offset, *it});
            }
        }
    }

    //Find cuts on each picked long read
    const int MIN_HITS = 0.25 * query_minimizers.size();
    vector<pair<int, cut> > cuts_tmp;
    for (auto it = hits.begin(); it != hits.end(); it++) {
        if (it->second.size() < MIN_HITS) {
            continue;
        }
        cut ans = find_range(it->second);
        if (ans.range.first == 0 && ans.range.second == 0)
            continue;
        ans.number_of_minimizers = it->second.size();
        cuts_tmp.push_back({it->first, ans});
    }

    vector<pair<string, pair<pair<int, int>, int> > > cuts;
    if (cuts_tmp.empty())
        cuts =  vector<pair<string, pair<pair<int, int>, int> > >();
    else if (classify)
        cuts = classify_reads(hits, cuts_tmp);
    else {
        for (int i = 0; i < cuts_tmp.size(); i++) {
            pair<int, cut> curr = cuts_tmp[i];
            cuts.push_back({sequences[curr.first].first, {curr.second.range, sequences[curr.first].second}});
        }
    }

    return cuts;
}