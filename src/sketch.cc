#include <deque>
#include <mutex>
#include <thread>
#include <string>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include "logger.h"
#include "common.h"
#include "MurmurHash3.h"

#include "sketch.h"

using namespace std;

//TODO FIX
const int MAXN = 30000;
const int BUCKET_SIZE = 100;

ProgressBar sketch_progress(80);
ProgressBar sort_progress(80);
uint64_t size_read = 0;
int sort_completed = 0;

mutex read_mutex, name_mutex, sort_mutex;

Sketch::Sketch() {}

std::ifstream::pos_type filesize(const string filename) {
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

Sketch::Sketch(string lp, string dp, int k, int w) {
    kmer_size = k;
    window_size = w;
    lr_path = lp;
	dat_path = dp;

    gz_fin = gzopen(lr_path.c_str(), "rb");
    if (!gz_fin) {
        Logger::instance().info("Could not create file. Quitting.\n");
        exit(1);
    }

    file_size = filesize(lr_path);
    zbuffer = (char *) malloc(BUFFSIZE);

	build_sketch_mt();
}

Sketch::Sketch(string dp) {
	dat_path = dp;
	load();
	compute_freq_th();
}

void Sketch::dump(vector<pair<uint64_t, Location> > &ref_minimizers_vec) {
    //Write minimizers
    ofstream fout(dat_path + "/sketch.dat", ios::out | ios::binary);
    if (!fout) {
        Logger::instance().info("Could not create file. Quitting.\n");
        exit(1);
    }

    ProgressBar write_prog(80);
    write_prog.update(0.0, "Writing Sketch ");

    uint32_t hash_entries = 1;
    uint64_t prv_hash, cnt = 0;
    vector<pair<uint64_t, uint32_t> > hash_sizes;
    uint64_t total_entries = ref_minimizers_vec.size();
    fout.write((char*)&total_entries, sizeof(uint64_t));

    int write_buf_size = 0;
    Location* write_buf = (Location*)malloc(1000000 * sizeof(Location));

    if (!ref_minimizers_vec.empty()) {
        prv_hash = ref_minimizers_vec[0].first;
        fout.write((char*)&(ref_minimizers_vec[0].second), sizeof(Location));

        for (int i = 1; i < ref_minimizers_vec.size(); i++) {
            if (ref_minimizers_vec[i].first == prv_hash)
                hash_entries++;
            else {
                hash_sizes.push_back({prv_hash, hash_entries});
                cnt += hash_entries;
                prv_hash = ref_minimizers_vec[i].first;
                hash_entries = 1;
            }
            if (write_buf_size == 1000000) {
                fout.write((char*)write_buf, 1000000 * sizeof(Location));
                write_buf_size = 0;
                write_prog.update(((float)cnt/(float)total_entries) * 100, "Writing Sketch ");
            }
            write_buf[write_buf_size++] = ref_minimizers_vec[i].second;
        }
        if (write_buf_size > 0)
            fout.write((char*)write_buf, write_buf_size * sizeof(Location));
        cnt += hash_entries;
        hash_sizes.push_back({prv_hash, hash_entries});

        for (int i = 0; i < hash_sizes.size(); i++) {
            fout.write((char*)&(hash_sizes[i].first), sizeof(uint64_t));
            fout.write((char*)&(hash_sizes[i].second), sizeof(uint32_t));
        }
    }

    write_prog.update(((float)cnt/(float)total_entries) * 100, "Writing Sketch ");

    fout.close();

    //Write sequences
    ofstream seqout(dat_path + "/sequences.dat", ios::out | ios::binary);
    for (uint32_t i = 0; i < names.size(); i++) {
        uint8_t size = names[i].first.length();
        seqout.write((char*)&i, sizeof(uint32_t));
        seqout.write((char*)&size, sizeof(uint8_t));
        seqout.write(names[i].first.c_str(), size);
        seqout.write((char*)&names[i].second, sizeof(uint16_t));
    }
    seqout.close();
}

void Sketch::load() {
    ifstream fin(dat_path + "/sketch.dat", ios::in | ios::binary);
	if (!fin) {
		Logger::instance().info("Could not load file. Quitting.\n");
		exit(1);
	}

	uint64_t cnt = 0;
	ProgressBar progress(80);
	progress.update(0.0, "Loading Sketch");

	uint64_t total_entries;
	fin.read(reinterpret_cast<char*>(&total_entries), sizeof(uint64_t));
    ref_minimizers.resize(total_entries);
    fin.read(reinterpret_cast<char*>(&ref_minimizers[0]), total_entries * sizeof(Location));

    uint64_t hash, offset = 0;
    uint32_t entries;
    while(fin.read(reinterpret_cast<char*>(&hash), sizeof(uint64_t))) {
        fin.read(reinterpret_cast<char*>(&entries), sizeof(uint32_t));
		cnt += entries;
        hashes.insert({hash, {offset, entries}});
        offset += entries;
		progress.update(((float)cnt / (float)total_entries) * 100, "Loading Sketch");
    }
    fin.close();

//    for (auto it = hashes.begin(); it != hashes.end(); it++) {
//        cout << it->first << ": ";
//        for (int i = it->second.first; i < it->second.first + it->second.second; i++)
//            cout << ref_minimizers[i].seq_id << "-" << ref_minimizers[i].offset << ", ";
//        cout << endl;
//    }

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
        names.push_back({name, len});
    }
    seqin.close();
}

void Sketch::sketch_query(vector<string> reads, int k, int w) {
    query_minimizers.clear();
    rev_query_minimizers.clear();
	kmer_size = k;
	window_size = w;
	build_sketch_query(reads);
}

void Sketch::sketch_query(string query, int k, int w) {
    query_minimizers.clear();
    rev_query_minimizers.clear();
    kmer_size = k;
    window_size = w;
    vector<string> queries;
    queries.push_back(query);
    build_sketch_query(queries);
}

void Sketch::read_buffer() {
    buff_size = gzread(gz_fin, zbuffer, BUFFSIZE);
    size_read = gzoffset(gz_fin);
    buff_pos = 0;
    if (buff_size != 0)
        sketch_progress.update(((float)size_read/(float)file_size) * 100, "Building Sketch");
    if (size_read == file_size)
        sort_progress.update(0.0, "Sorting Entries");
    if (buff_size == 0 and gzeof(gz_fin) == 0) {
        buff_size = -1;
    }
    if (buff_size < 0) {
        int err;
        fprintf(stderr, "gzread error: %s\n", gzerror(gz_fin, &err));
        exit(1);
    }
}

inline uint32_t Sketch::read_line(string& seq) {
    char cur;

    uint32_t i = 0;
    while (true) {
        if (buff_pos >= buff_size) {
            read_buffer();
            if (buff_size == 0)
                return 0;
        }

        cur = zbuffer[buff_pos++];
        if (cur == '\n') {
            return i;
        }

        seq += cur;
        i++;
    }
}

void Sketch::build_sketch(int tid, ProgressBar progress, vector<pair<uint64_t, Location> > &ref_minimizers_vec) {
    string name, read;
    name.reserve(100);
    read.reserve(100000);
    int curr_id;
    uint16_t nm_sz;
    int rd_sz;

    vector<pair<uint64_t, Location> > new_vec;
    new_vec.reserve(10000000);
    while (true) {
        // MUTEX
        read_mutex.lock();
        nm_sz = read_line(name);
        if (!nm_sz) {
            read_mutex.unlock();
            break;
        }
        rd_sz = read_line(read);
        curr_id = read_id;
        read_id++;
        read_mutex.unlock();
        // MUTEX END

        transform(read.begin(), read.end(), read.begin(), ::toupper);
        get_ref_minimizers(&read[0], curr_id, read.size(), new_vec);

        name_mutex.lock();
        names.push_back({name.substr(1, nm_sz - 1), rd_sz});
        name_mutex.unlock();

        read.clear();
        name.clear();
    }

    sort(new_vec.begin(), new_vec.end(), [](const auto &a, const auto &b) { return a.first < b.first;});

    sort_mutex.lock();
    sort_completed++;
    sort_progress.update(((float)sort_completed/(float)thread_cnt) * 100, "Sorting Entries");
    sort_mutex.unlock();

    ref_minimizers_vec = new_vec;
}

void Sketch::build_sketch_mt() {
	sketch_progress.update(0.0, "Building Sketch");

    vector<pair<uint64_t, Location> > ref_minimizers_vec[thread_cnt];

	thread myThreads[thread_cnt];
    for (int i = 0; i < thread_cnt; i++){
        myThreads[i] = thread(&Sketch::build_sketch, this, i, sketch_progress, std::ref(ref_minimizers_vec[i]));
    }
    for (int i = 0; i < thread_cnt; i++){
        myThreads[i].join();
    }

    thread merge_threads[thread_cnt/2];
    int mx_th = thread_cnt/2;
    ProgressBar merge_prog(80);
    merge_prog.update(0.0, "Merging Entries");
    while (mx_th > 0) {
        for (int i = 0; i < mx_th; i++) {
            merge_threads[i] = thread(&Sketch::merge, this, ref(ref_minimizers_vec[i]), ref(ref_minimizers_vec[mx_th*2 - i - 1]));
        }
        for (int i = 0; i < mx_th; i++){
            merge_threads[i].join();
        }
        merge_prog.update((1.0/(float)mx_th) * 100, "Merging Entries");
        mx_th /= 2;
    }

    dump(ref_minimizers_vec[0]);
}

void Sketch::merge(vector<pair<uint64_t, Location> > &a, vector<pair<uint64_t, Location> > &b) {
    int mid = a.size();
    a.insert(a.end(), b.begin(), b.end());
    b.clear();
    b.resize(0);
    inplace_merge(a.begin(), a.begin() + mid, a.end(),
                  [](auto &a, auto &b) {return a.first < b.first;});
}

void Sketch::build_sketch_query(vector<string> reads) {
    for (int i = 0; i < reads.size(); i++) {
        transform(reads[i].begin(), reads[i].end(), reads[i].begin(), ::toupper);
        get_query_minimizers(&reads[i][0], i, reads[i].size());
    }
    update_query_sketch(FRW);
    query_minimizers_vec.clear();

    for (int i = 0; i < reads.size(); i++) {
        string q = reverse_complement(reads[i]);
        get_query_minimizers(&q[0], i, q.size());
    }
    update_query_sketch(REV);
    query_minimizers_vec.clear();
}

void Sketch::get_ref_minimizers(char* read, int id, int len, vector<pair<uint64_t, Location> > &ref_minimizers_vec) {
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
            if (ref_minimizers_vec.empty() || ref_minimizers_vec.back().first != window.front().first) {
                ref_minimizers_vec.push_back(window.front());
            }
        }
    }
}

void Sketch::get_query_minimizers(char* read, int id, int len) {
    deque<pair<uint64_t, int> > window;

    char tmp[16];

    int win = window_size - 1;
    for (int i = 0; i < len - kmer_size + 1; i++) {

        MurmurHash3_x64_128(read + i, kmer_size, 50, tmp);
        uint64_t mhash = *((uint64_t *)tmp);

        while(!window.empty() && window.back().first >= mhash)
            window.pop_back();

        window.push_back(pair<uint64_t, int>(mhash, i));

        while (!window.empty() && window.front().second <= i - window_size)
            window.pop_front();

        if (i >= win) {
            if (query_minimizers_vec.empty() || query_minimizers_vec.back().first != window.front().first)
                query_minimizers_vec.push_back(window.front());
        }
    }
}

void Sketch::update_query_sketch(int ort) {
    for (int i = 0; i < query_minimizers_vec.size(); i++) {
        if (ort == FRW)
            query_minimizers.insert(query_minimizers_vec[i].first);
        else
            rev_query_minimizers.insert(query_minimizers_vec[i].first);
    }
}

void Sketch::compute_freq_th() {
    vector<pair<uint64_t, int> > frequencies;
    for (auto it = hashes.begin(); it != hashes.end(); it++) {
        int tmp_frq = it->second.second;
        if (tmp_frq > freq_th) {
            freq_th = tmp_frq;
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

vector<pair<string, pair<pair<int, int>, pair<int, int> > > > Sketch::classify_reads(map<int, vector<pair<int, uint64_t> > > hits,
                                                                         vector<pair<int, cut> > cuts_tmp) {
    if (cuts_tmp.empty())
        return vector<pair<string, pair<pair<int, int>, pair<int, int> > > >();

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
    vector<pair<string, pair<pair<int, int>, pair<int, int> > > > cuts;

    int i = 0;
    //Separate left and right minimizers
    while (i < cuts_tmp.size() && cuts_tmp[i].second.type == BIMODAL) {
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
        cuts.push_back({names[curr_read].first, {curr_cut.range, {BIMODAL, curr_cut.orientation}}});
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
                cuts.push_back({names[curr_read].first, {curr_cut.range, {LEFT, curr_cut.orientation}}});
            else if (right_similarity > left_similarity)
                cuts.push_back({names[curr_read].first, {curr_cut.range, {RIGHT, curr_cut.orientation}}});
            else
                cuts.push_back({names[curr_read].first, {curr_cut.range, {MISC, curr_cut.orientation}}});
        }
        //Case 2: No bimodal reads, must classify reads on the go -> first set of seen minimizers are assigned to left
        else if (left_minimizers.empty()) {
            for (int j = 0; j < hits[curr_read].size(); j++) {
                left_minimizers.insert(hits[curr_read][j].second);
            }
            cuts.push_back({names[curr_read].first, {curr_cut.range, {LEFT, curr_cut.orientation}}});
        }
        //Case 3: all seen reads have been classified as left and right is still empty -> if this new read doesn't have
        //        enough similarity with the left set, assign it to right
        else if (right_minimizers.empty()) {
            float left_similarity = minimizer_similarity(left_minimizers, hits[curr_read]);
            if (left_similarity > 0.5) {
                cuts.push_back({names[curr_read].first, {curr_cut.range, {LEFT, curr_cut.orientation}}});
                continue;
            }
            for (int j = 0; j < hits[curr_read].size(); j++) {
                right_minimizers.insert(hits[curr_read][j].second);
            }
            cuts.push_back({names[curr_read].first, {curr_cut.range, {RIGHT, curr_cut.orientation}}});
        }
    }
    return cuts;
}

vector<pair<string, pair<pair<int, int>, pair<int, int> > > > Sketch::find_cuts(bool classify) {
    map<int, vector<pair<int, uint64_t> > > frw_hits;
    map<int, vector<pair<int, uint64_t> > > rev_hits;

    vector<pair<string, pair<pair<int, int>, pair<int, int> > > > cuts, cuts_frv, cuts_rev;

    //Find all forward hits
    for (auto it = query_minimizers.begin(); it != query_minimizers.end(); it++) {
        auto idx = hashes.find(*it);
        if (idx == hashes.end())
            continue;
        if (idx->second.second >= freq_th)
            continue;
        for (int i = idx->second.first; i < idx->second.first + idx->second.second; i++) {
            auto j = frw_hits.find(ref_minimizers[i].seq_id);
            if (j == frw_hits.end()) {
                vector<pair<int, uint64_t> > tmp;
                tmp.reserve(32);
                tmp.push_back({ref_minimizers[i].offset, *it});
                frw_hits.insert({ref_minimizers[i].seq_id, tmp});
            }
            else {
                j->second.push_back({ref_minimizers[i].offset, *it});
            }
        }
    }

    for (auto it = rev_query_minimizers.begin(); it != rev_query_minimizers.end(); it++) {
        auto idx = hashes.find(*it);
        if (idx == hashes.end())
            continue;
        if (idx->second.second >= freq_th)
            continue;
        for (int i = idx->second.first; i < idx->second.first + idx->second.second; i++) {
            auto j = rev_hits.find(ref_minimizers[i].seq_id);
            if (j == rev_hits.end()) {
                vector<pair<int, uint64_t> > tmp;
                tmp.reserve(32);
                tmp.push_back({ref_minimizers[i].offset, *it});
                rev_hits.insert({ref_minimizers[i].seq_id, tmp});
            }
            else {
                j->second.push_back({ref_minimizers[i].offset, *it});
            }
        }
    }

    //Find cuts on each picked long read
    int MIN_HITS = 0.25 * query_minimizers.size();;

    vector<pair<int, cut> > cuts_tmp_frw;
    for (auto it = frw_hits.begin(); it != frw_hits.end(); it++) {
        auto r = rev_hits.find(it->first);
        if (r != rev_hits.end())
            continue;
        if (it->second.size() < MIN_HITS) {
            continue;
        }
        cut ans = find_range(it->second);
        ans.orientation = FRW;
        if (ans.range.first == 0 && ans.range.second == 0)
            continue;
        ans.number_of_minimizers = it->second.size();
        cuts_tmp_frw.push_back({it->first, ans});
    }
	if (cuts_tmp_frw.size() > 200)
                return cuts;

    MIN_HITS = 0.25 * rev_query_minimizers.size();
    vector<pair<int, cut> > cuts_tmp_rev;
    for (auto it = rev_hits.begin(); it != rev_hits.end(); it++) {
        auto f = frw_hits.find(it->first);
        if (f != frw_hits.end())
            continue;
        if (it->second.size() < MIN_HITS) {
            continue;
        }
        cut ans = find_range(it->second);
        ans.orientation = REV;
        if (ans.range.first == 0 && ans.range.second == 0)
            continue;
        ans.number_of_minimizers = it->second.size();
        cuts_tmp_rev.push_back({it->first, ans});
    }
	if (cuts_tmp_frw.size() + cuts_tmp_rev.size() > 200)
                return cuts;

    if (cuts_tmp_frw.empty() && cuts_tmp_rev.empty())
        cuts =  vector<pair<string, pair<pair<int, int>, pair<int, int> > > >();
    else if (classify) {
        cuts_frv = classify_reads(frw_hits, cuts_tmp_frw);
        cuts_rev = classify_reads(rev_hits, cuts_tmp_rev);
        for (int i = 0; i < cuts_frv.size(); i++) {
            cuts.push_back(cuts_frv[i]);
        }
        for (int i = 0; i < cuts_rev.size(); i++) {
            cuts.push_back(cuts_rev[i]);
        }
    }
    else {
        for (int i = 0; i < cuts_tmp_frw.size(); i++) {
            pair<int, cut> curr = cuts_tmp_frw[i];
            cuts.push_back({names[curr.first].first, {curr.second.range, {names[curr.first].second, curr.second.orientation}}});
        }
    }

    return cuts;
}

void fix_reverse(vector<pair<pair<string, string>, pair<pair<int, int>, pair<int, int> > > > &cuts) {
    aligner al(30000);

    vector<int> frw_bimodal, rev_bimodal, frw_right, rev_right, frw_left, rev_left;

    for (int i = 0; i < cuts.size(); i++) {
        if (cuts[i].second.second.second == FRW) {
            if (cuts[i].second.second.first == BIMODAL)
                frw_bimodal.push_back(i);
            else if (cuts[i].second.second.first == LEFT)
                frw_left.push_back(i);
            else if (cuts[i].second.second.first == RIGHT)
                frw_right.push_back(i);
        }
        if (cuts[i].second.second.second == REV) {
            if (cuts[i].second.second.first == BIMODAL)
                rev_bimodal.push_back(i);
            else if (cuts[i].second.second.first == LEFT)
                rev_left.push_back(i);
            else if (cuts[i].second.second.first == RIGHT)
                rev_right.push_back(i);
        }
    }

    if ((rev_bimodal.empty() && rev_left.empty() && rev_left.empty()) || (frw_bimodal.empty() && frw_left.empty() && frw_right.empty()))
        return;

    if (!rev_bimodal.empty()) {
        for (int i = 0; i < rev_left.size(); i++)
            cuts[rev_left[i]].second.second.first = RIGHT;
        for (int i = 0; i < rev_right.size(); i++)
            cuts[rev_right[i]].second.second.first = LEFT;
    }

    //Case 1: rev and frw have bimodal reads
    if (!frw_bimodal.empty() && !rev_bimodal.empty())
        return;

    //Case 2: frw has bimodal, rev does not
    else if (!frw_bimodal.empty() && rev_bimodal.empty()) {
        string rev_1 = cuts[rev_left[0]].first.second;
        if (!frw_left.empty()) {
            string left = cuts[frw_left[0]].first.second;

            al.align(left, rev_1);
            if (al.get_left_anchor() < 20 && al.get_right_anchor() < 20) {
                for (int i = 0; i < rev_left.size(); i++)
                    cuts[rev_left[i]].second.second.first = RIGHT;
                for (int i = 0; i < rev_right.size(); i++)
                    cuts[rev_right[i]].second.second.first = LEFT;
            }
        }
    }

    //Case 3: rev has bimodal, frw does not
    else if (frw_bimodal.empty() && !rev_bimodal.empty()) {
        if (!rev_left.empty()) {
            string right = cuts[rev_left[0]].first.second;
            string frw_1 = cuts[frw_left[0]].first.second;

            al.align(right, frw_1);
            if (al.get_left_anchor() > 20 || al.get_right_anchor() > 20) {
                for (int i = 0; i < frw_left.size(); i++)
                    cuts[frw_left[i]].second.second.first = RIGHT;
                for (int i = 0; i < frw_right.size(); i++)
                    cuts[frw_right[i]].second.second.first = LEFT;
            }
        }
    }

    else {
        //Case 4: single peak
        if (frw_right.empty() && rev_right.empty())
            return;

        //Case 5: long insertion
        else {
            string frw_1 = cuts[frw_left[0]].first.second;
            string rev_1 = cuts[rev_left[0]].first.second;

            al.align(frw_1, rev_1);
            if (al.get_left_anchor() < 20 && al.get_right_anchor() < 20) {
                for (int i = 0; i < rev_left.size(); i++)
                    cuts[rev_left[i]].second.second.first = RIGHT;
                for (int i = 0; i < rev_right.size(); i++)
                    cuts[rev_right[i]].second.second.first = LEFT;
            }
        }
    }
}