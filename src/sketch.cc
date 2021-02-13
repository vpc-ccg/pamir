#include <cmath>
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

#include <cstring>

using namespace std;

//TODO FIX
const int MAXN = 30000;
const int BUCKET_SIZE = 100;
int REF_MIN = 10;

ProgressBar sketch_progress(80);
ProgressBar sort_progress(80);
uint64_t size_read = 0;
int sort_completed = 0;

mutex read_mutex, sort_mutex;

Sketch::Sketch() {}

ifstream::pos_type filesize(const string filename) {
    ifstream in(filename, ifstream::ate | ifstream::binary);
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

    build_sketch();
}

Sketch::Sketch(string dp, int len, int distance, int k, int w) {
    short_read_len = len;
    genome_anchor_distance = distance;
    kmer_size = k;
    window_size = w;
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

    uint64_t hash_entries = 1;
    hash_t prv_hash;
    vector<minimizer> hash_sizes;
    uint64_t total_entries = ref_minimizers_vec.size(), cnt = 0;
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
                hash_sizes.push_back((minimizer){.hash_value = prv_hash, .offset = cnt});
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

        hash_sizes.push_back((minimizer){.hash_value = prv_hash, .offset = cnt});
        cnt += hash_entries;

        uint64_t total_hashes = hash_sizes.size();
        fout.write((char*)&total_hashes, sizeof(uint64_t));
        fout.write((char*)&hash_sizes[0], total_hashes * sizeof(minimizer));
    }

    write_prog.update(((float)cnt/(float)total_entries) * 100, "Writing Sketch ");

    fout.close();

    //Write sequences
    ofstream seqout(dat_path + "/sequences.dat", ios::out | ios::binary);
    for (uint32_t i = 0; i < sequences.size(); i++) {
        uint8_t size = sequences[i].first.size();
        seqout.write((char*)&size, sizeof(uint8_t));
        seqout.write(&sequences[i].first[0], size);
        seqout.write((char*)&sequences[i].second, sizeof(offset_t));
    }
    seqout.close();
}

void Sketch::load() {
    ifstream fin(dat_path + "/sketch.dat", ios::in | ios::binary);
    if (!fin) {
        Logger::instance().info("Could not load file. Quitting.\n");
        exit(1);
    }

    uint64_t total_entries, total_hashes;

    //Read all locations
    fin.read(reinterpret_cast<char*>(&total_entries), sizeof(uint64_t));
    ref_minimizers.resize(total_entries);
    fin.read(reinterpret_cast<char*>(&ref_minimizers[0]), total_entries * sizeof(Location));

    //Read all hashes and sizes
    fin.read(reinterpret_cast<char*>(&total_hashes), sizeof(uint64_t));
    minimizers.resize(total_hashes);
    fin.read(reinterpret_cast<char*>(&minimizers[0]), total_hashes * sizeof(minimizer));

    fin.close();

    string seq_file = dat_path + "/sequences.dat";
    Logger::instance().info("Loading sequence names from: %s\n", seq_file.c_str());
    ifstream seqin(seq_file, ios::in | ios::binary);
    if (!seqin) {
        Logger::instance().info("Could not load file. Quitting.\n");
        exit(1);
    }

    string name;
    offset_t len;
    uint8_t name_size;
    while(seqin.read(reinterpret_cast<char*>(&name_size), sizeof(uint8_t))) {
        name.resize(name_size);
        seqin.read(&name[0], name_size);
        seqin.read(reinterpret_cast<char*>(&len), sizeof(offset_t));
        sequences.push_back({name, len});
    }
    seqin.close();
}

vector<cut> Sketch::query(vector<string> &reads, bool classify) {
    unordered_set<hash_t> frw_minimizers;
    unordered_set<hash_t> rev_minimizers;
    vector<pair<uint64_t, int> > minimizers_vec;

    build_query_sketch(reads, minimizers_vec, frw_minimizers, rev_minimizers);

    return find_cuts(classify, frw_minimizers, rev_minimizers);
}

vector<cut> Sketch::query(vector<string> &reads, bool classify, string& ref_l, string& ref_r) {
    unordered_set<hash_t> frw_minimizers;
    unordered_set<hash_t> rev_minimizers;
    vector<pair<uint64_t, int> > minimizers_vec;

    build_query_sketch(reads, minimizers_vec, frw_minimizers, rev_minimizers);

    return find_cuts(classify, frw_minimizers, rev_minimizers, ref_l, ref_r);
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

void Sketch::build_sketch_mt(int id, const ProgressBar progress, vector<pair<hash_t, Location> > &ref_minimizers_vec) {
    string name, read;
    name.reserve(100);
    read.reserve(100000);
    int curr_id;

    vector<pair<hash_t, Location> > new_vec;
    new_vec.reserve(10000000);
    while (true) {
        // MUTEX
        read_mutex.lock();
        name.clear();
        read.clear();
        read_line(name);
        if (!name.size()) {
            read_mutex.unlock();
            break;
        }
        read_line(read);
        sequences.push_back({name.substr(1), read.size()});
        curr_id = read_id;
        read_id++;
        read_mutex.unlock();
        // MUTEX END

        transform(read.begin(), read.end(), read.begin(), ::toupper);
        get_ref_minimizers(&read[0], curr_id, read.size(), new_vec);

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

void Sketch::build_sketch() {
    sketch_progress.update(0.0, "Building Sketch");

    vector<pair<hash_t, Location> > ref_minimizers_vec[thread_cnt];

    thread myThreads[thread_cnt];
    for (int i = 0; i < thread_cnt; i++){
        myThreads[i] = thread(&Sketch::build_sketch_mt, this, i, sketch_progress, ref(ref_minimizers_vec[i]));
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

void Sketch::merge(vector<pair<hash_t, Location> > &a, vector<pair<hash_t, Location> > &b) {
    int mid = a.size();
    a.insert(a.end(), b.begin(), b.end());
    b.clear();
    b.resize(0);
    inplace_merge(a.begin(), a.begin() + mid, a.end(),
                  [](auto &a, auto &b) {return a.first < b.first;});
}

void Sketch::build_query_sketch(vector<string> &reads, vector<pair<uint64_t, int> > &minimizers_vec,
                                unordered_set<hash_t> &frw_minimizers, unordered_set<hash_t> &rev_minimizers) {
    for (int i = 0; i < reads.size(); i++) {
        transform(reads[i].begin(), reads[i].end(), reads[i].begin(), ::toupper);
        get_query_minimizers(&reads[i][0], i, reads[i].size(), minimizers_vec);
    }
    for (int i = 0; i < minimizers_vec.size(); i++)
        frw_minimizers.insert(minimizers_vec[i].first);
    minimizers_vec.clear();

    for (int i = 0; i < reads.size(); i++) {
        string q = reverse_complement(reads[i]);
        get_query_minimizers(&q[0], i, q.size(), minimizers_vec);
    }
    for (int i = 0; i < minimizers_vec.size(); i++)
        rev_minimizers.insert(minimizers_vec[i].first);
    minimizers_vec.clear();
}

void Sketch::get_ref_minimizers(char* read, id_t id, int len, vector<pair<hash_t, Location> > &ref_minimizers_vec) {
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

void Sketch::get_query_minimizers(char* read, id_t id, offset_t len, vector<pair<uint64_t, int> > &minimizers_vec) {
    deque<pair<uint64_t, int> > window;

    char tmp[16];

    int win = window_size - 1;
    for (auto i = 0; i < len - kmer_size + 1; i++) {

        MurmurHash3_x64_128(read + i, kmer_size, 50, tmp);
        uint64_t mhash = *((uint64_t *)tmp);

        while(!window.empty() && window.back().first >= mhash)
            window.pop_back();

        window.push_back(pair<uint64_t, int>(mhash, i));

        while (!window.empty() && window.front().second <= i - window_size)
            window.pop_front();

        if (i >= win) {
            if (minimizers_vec.empty() || minimizers_vec.back().first != window.front().first) {
                minimizers_vec.push_back(window.front());
            }
        }
    }
}

void Sketch::compute_freq_th() {
    vector<hash_size_t> frequencies;
    hash_t prv_offset = 0;
    for (auto it = minimizers.begin() + 1; it != minimizers.end(); it++) {
        hash_size_t frq = it->offset - prv_offset;
        prv_offset = it->offset;
        frequencies.push_back(frq);
    }
    frequencies.push_back(ref_minimizers.size() - prv_offset);
    sort(frequencies.begin(), frequencies.end(), [](auto &a, auto &b) {
        return a > b;
    });
    hash_size_t top = minimizers.size() * 0.001;
    freq_th = frequencies[top];
}

int get_avg(vector<hash_size_t> v) {
    hash_size_t sum = 0;
    for (int i = 0; i < v.size(); i++) {
        sum += v[i];
    }
    int res = sum / v.size();
    return res;
}

int get_avg_new(vector<float> v) {
    float sum = 0;
    for (int i = 0; i < v.size(); i++) {
        sum += v[i];
    }
    int res = sum / v.size();
    return res;
}

cut Sketch::find_range(vector<hit>& hits, mem_offset_t start, hash_size_t size) {
    sort(hits.begin() + start, hits.begin() + start + size, hit());

    //TODO: Optimize

    offset_t min_offset = hits[start].offset;
    offset_t max_offset = hits[start + size - 1].offset + kmer_size;
    int n = ceil(((max_offset - min_offset)/100) + 2) * 100;

    //Create buckets
    vector<pair<offset_t, offset_t> > buckets;
    buckets.reserve(BUCKET_SIZE);
    map<int, float> cnt_new;
    for (int s = max(min_offset - 200, 0); s < min_offset + n; s += BUCKET_SIZE) {
        buckets.push_back(make_pair(s, s + BUCKET_SIZE));
    }

    vector<pair<offset_t, offset_t> > locs;
    for (int i = 0; i < buckets.size(); i++) {
        locs.push_back({buckets[i].second, buckets[i].first});
    }

     //Populate buckets by hit locations
    cnt_new[0] = 0;
    cnt_new[buckets.size() - 1] = 0;
    int j = 1;
    float weight = 0;
    for (auto i = start; i < start + size; i++) {
        offset_t curr_hit_start = hits[i].offset;
        offset_t curr_hit_end = hits[i].offset + kmer_size;
        float l_weight, r_weight;
        while (curr_hit_start >= buckets[j].second) {
            j++;
        }
        if (curr_hit_end >= buckets[j].second) {
            l_weight = (buckets[j].second - curr_hit_start)/(float)kmer_size;
            r_weight = 1 - l_weight;
            locs[j].second = buckets[j].second;
            locs[j+1].second = max((int)locs[j+1].second, (int)curr_hit_end);
            locs[j+1].first = buckets[j+1].first;
        }
        else {
            l_weight = 1.0;
            r_weight = 0;
            locs[j].second = max((int)locs[j].second, (int)curr_hit_end);
        }
        cnt_new[j] += l_weight;
        cnt_new[j + 1] += r_weight;
        locs[j].first = min((int)locs[j].first, (int)curr_hit_start);
    }

    //Find bucket weights
    vector<float> freq_new;
    freq_new.reserve(BUCKET_SIZE);
    for (int i = 0; i < buckets.size(); i++) {
#ifdef DEBUG
        Logger::instance().debug("%d-%d: %f\n", buckets[i].first, buckets[i].second, cnt_new[i]);
#endif
        freq_new.push_back(cnt_new[i]);
    }
    //find cut around highest peak
    int mx_idx = max_element(freq_new.begin(), freq_new.end()) - freq_new.begin();

    int avg = get_avg_new(freq_new);

    offset_t l = mx_idx;
    while (l > 0) {
        if (freq_new[l - 1] > 0) {
            l--;
        } else {
            break;
        }
    }
    offset_t r = mx_idx;
    while (r < (int)freq_new.size() && freq_new[r] > 0) {
        r++;
    }
    range_s range_1 = {l, r};

    //find cut around second highest peak
    range_s range_2;
    //find if peak is on the left side or the right side of maximum peak

    float r_w, l_w;

    //TODO: check sum of weights. not just the one bucket
    int mx_idx_l = max_element(freq_new.begin(), freq_new.begin() + max(l - 1, 0)) - freq_new.begin();
    int mx_idx_r = max_element(freq_new.begin() + min(r + 1, (int) freq_new.size() - 1), freq_new.end()) - freq_new.begin();

    l_w = freq_new[mx_idx_l];
    r_w = freq_new[mx_idx_r];

    if (l - 1 <= 0) {
        l_w = -1;
    }
    if (r + 1 >= freq_new.size()) {
        r_w = -1;
    }

    if (l_w == r_w) {
        range_2 = {0, 0};
    }
    else {
        //    float prv_peak = freq_new[mx_idx] + freq_new[mx_idx - 1] + freq_new[mx_idx + 1];
        float prv_peak = freq_new[mx_idx];
        mx_idx = l_w > r_w ? mx_idx_l : mx_idx_r;
//    float new_peak = freq_new[mx_idx] + freq_new[mx_idx - 1] + freq_new[mx_idx + 1];
        float new_peak = freq_new[mx_idx];

        if (mx_idx <= r && mx_idx_l >= l)
            range_2 = {0, 0};
            //Discard if new peak is too low comparing to the maximum peak (probably doesn't belong to this insertion)
        else if (new_peak < prv_peak/2) {
            range_2 = {0, 0};
        }
        else {
            l = mx_idx;
            while (l > 0) {
                if (freq_new[l - 1] > 0) {
                    l--;
                } else {
                    break;
                }
            }
            r = mx_idx;
            while (r < (int)freq_new.size() && freq_new[r] > 0) {
                r++;
            }
            range_2 = {l, r};
        }
    }

    range_s final_range;
    uint8_t type;
    cut ans;

    if (range_2.start == range_2.end) {
        final_range.start = max(0, locs[range_1.start].first - 300);
        final_range.end = locs[range_1.end - 1].second + 300;
        ans.type = SINGLE_PEAK;
        ans.peak1 = final_range;
        ans.peak2 = {0, 0};
    }
        //Case 2: Two peaks found/Bimodal read
    else {
        ans.type = BIMODAL;
        //If second peak is on the right side of the first peak
        if (range_1.start < range_2.start) {
            final_range.start = max(0, locs[range_1.start].first - 300);
            final_range.end = locs[range_2.end - 1].second + 300;
            ans.peak1 = {locs[range_1.start].first, locs[range_1.end - 1].second};
            ans.peak2 = {locs[range_2.start].first, locs[range_2.end - 1].second};
        }
        //If second peak is on the left side of the first peak
        else {
            final_range.start = max(0, locs[range_2.start].first - 300);
            final_range.end = locs[range_1.end - 1].second + 300;
            ans.peak2 = {locs[range_1.start].first, locs[range_1.end - 1].second};
            ans.peak1 = {locs[range_2.start].first, locs[range_2.end - 1].second};
        }
        ans.estimated_insertion = abs(ans.peak1.start - ans.peak2.start);
    }
    ans.range = final_range;

    return ans;
}

pair<float, float> Sketch::minimizer_similarity(unordered_set<hash_t>& ref, vector<hit>& q, int start, int size) {
    int cnt = 0;
    for (auto i = start; i < start + size; i++) {
        if (find(ref.begin(), ref.end(), q[i].hash_value) != ref.end())
            cnt++;
    }
    return {float(cnt)/float(size), float(cnt)/float(ref.size())};
}

float Sketch::minimizer_similarity_single(unordered_set<hash_t>& ref, unordered_set<hash_t>& q) {
    int cnt = 0;
    for (auto i = q.begin(); i != q.end(); i++) {
        if (find(ref.begin(), ref.end(), *i) != ref.end())
            cnt++;
    }
    return float(cnt)/float(q.size());
}

float Sketch::compare_sequences(string& seq_a, string& seq_b) {
    unordered_set<hash_t> minimizers_a;
    unordered_set<hash_t> minimizers_b;
    vector<pair<uint64_t, int> > minimizers_vec;

    //TODO: move upper to extraction
    transform(seq_a.begin(), seq_a.end(), seq_a.begin(), ::toupper);
    get_query_minimizers(&seq_a[0], 0, seq_a.size(), minimizers_vec);

    for (int i = 0; i < minimizers_vec.size(); i++)
        minimizers_a.insert(minimizers_vec[i].first);
    minimizers_vec.clear();

    transform(seq_b.begin(), seq_b.end(), seq_b.begin(), ::toupper);
    get_query_minimizers(&seq_b[0], 0, seq_b.size(), minimizers_vec);

    for (int i = 0; i < minimizers_vec.size(); i++)
        minimizers_b.insert(minimizers_vec[i].first);
    minimizers_vec.clear();

    return minimizer_similarity_single(minimizers_a, minimizers_b);
}

void Sketch::classify_reads(vector<hit> &hits, vector<cut> &cuts) {
    if (cuts.empty())
        return;
    //Bimodal reads first, then ordered by number of minimizers used to pick the reads
    sort(cuts.begin(), cuts.end());

    unordered_set<hash_t> left_minimizers, right_minimizers;

    int i = 0;
    //Separate left and right minimizers
    while (i < cuts.size() && cuts[i].type == BIMODAL) {
        id_t curr_read = cuts[i].seq_id;

        auto h = equal_range(hits.begin(), hits.end(), curr_read, hit());
        for (auto it = h.first; it != h.second; it++) {
            //Left minimizers
            if (it->offset < cuts[i].peak1.end)
                left_minimizers.insert(it->hash_value);
                //Right minimizers
            else if (it->offset > cuts[i].peak2.start)
                right_minimizers.insert(it->hash_value);
        }
        i++;
    }

    for (i; i < cuts.size(); i++) {
        id_t curr_read = cuts[i].seq_id;

        auto h = equal_range(hits.begin(), hits.end(), curr_read, hit());

        int start = h.first - hits.begin();
        int size = h.second - h.first;
        //Case 1: We have seen a bimodal read or we have already populated the left and right sets
        // -> just compare similarity to left and right sets
        if (!left_minimizers.empty() && !right_minimizers.empty()) {
            pair<float, float> left_similarity = minimizer_similarity(left_minimizers, hits, start, size);
            pair<float, float> right_similarity = minimizer_similarity(right_minimizers, hits, start, size);

            if (abs(left_similarity.first - left_similarity.second) > 0.3 ||
                abs(right_similarity.first - right_similarity.second) > 0.3)
                cuts[i].type = MISC;
            else if (left_similarity > right_similarity) {
                cuts[i].type = PARTIAL_LEFT;
            }
            else if (right_similarity > left_similarity) {
                cuts[i].type = PARTIAL_RIGHT;
            }
            else {
                cuts[i].type = MISC;
            }
        }
            //Case 2: No bimodal reads, must classify reads on the go -> first set of seen minimizers are assigned to left
        else if (left_minimizers.empty()) {
            for (auto j = h.first; j != h.second; j++) {
                left_minimizers.insert(j->hash_value);
            }
            cuts[i].type = PARTIAL_LEFT;
        }
            //Case 3: all seen reads have been classified as left and right is still empty -> if this new read doesn't have
            //        enough similarity with the left set, assign it to right
        else if (right_minimizers.empty()) {
            pair<float, float> left_similarity = minimizer_similarity(left_minimizers, hits, start, size);
            if (abs(left_similarity.first - left_similarity.second) > 0.1) {
                cuts[i].type = MISC;
                continue;
            }
            if (left_similarity.first > 0.5) {
                cuts[i].type = PARTIAL_LEFT;
                continue;
            }
            for (auto j = h.first; j != h.second; j++) {
                right_minimizers.insert(j->hash_value);
            }
            cuts[i].type = PARTIAL_RIGHT;
        }
    }
    return;
}

pair<mem_offset_t, hash_size_t> Sketch::find_hit(const hash_t &hv) {
    vector<minimizer>::iterator it = lower_bound(minimizers.begin(), minimizers.end(), hv, minimizer());

    if (it == minimizers.end()) {
        return {0, 0};
    }
    else if (it->hash_value == hv) {
        mem_offset_t offset = it->offset;
        hash_size_t size = 0;
        if (it+1 != minimizers.end()) {
            size = (it+1)->offset - it->offset;
        }
        else {
            size = ref_minimizers.size() - it->offset;
        }
        return {offset, size};
    }
    else {
        return {0, 0};
    }
};

pair<vector<hit>, vector<hit> > Sketch::get_hits(unordered_set<hash_t>& query_frw, unordered_set<hash_t>& query_rev) {
    vector<hit> hits_frw;
    vector<hash_t> frequent_hashes;
    for (auto it = query_frw.begin(); it != query_frw.end(); it++) {
        auto idx = find_hit(*it);

        if (idx.second == 0 || idx.second >= freq_th) {
            frequent_hashes.push_back(*it);
            continue;
        }

        for (auto i = idx.first; i < idx.first + idx.second; i++) {
            hits_frw.push_back((hit) {.seq_id = ref_minimizers[i].seq_id, .offset = ref_minimizers[i].offset, .hash_value = *it});
        }
    }
    sort(hits_frw.begin(), hits_frw.end());
    for (int i = 0; i < frequent_hashes.size(); i++)
        query_frw.erase(frequent_hashes[i]);
    frequent_hashes.clear();


    vector<hit> hits_rev;
    for (auto it = query_rev.begin(); it != query_rev.end(); it++) {
        auto idx = find_hit(*it);

        if (idx.second == 0 || idx.second >= freq_th) {
            frequent_hashes.push_back(*it);
            continue;
        }

        for (auto i = idx.first; i < idx.first + idx.second; i++) {
            hits_rev.push_back((hit) {.seq_id = ref_minimizers[i].seq_id, .offset = ref_minimizers[i].offset, .hash_value = *it});
        }
    }
    sort(hits_rev.begin(), hits_rev.end());
    for (int i = 0; i < frequent_hashes.size(); i++)
        query_rev.erase(frequent_hashes[i]);
    frequent_hashes.clear();

    return {hits_frw, hits_rev};
}

pair<vector<hit>, vector<hit> > Sketch::get_hits(vector<pair<uint64_t, int> >& query_frw, vector<pair<uint64_t, int> >& query_rev) {
    vector<hit> hits_frw;
    vector<pair<hash_t, int> > hashes;

    hash_t prv_hash = UINT64_MAX;

    for (auto it = query_frw.begin(); it != query_frw.end(); it++) {
        auto idx = find_hit(it->first);

        if (idx.second == 0 || idx.second >= freq_th) {
            continue;
        }

        hashes.push_back(*it);

        if (it->first != prv_hash) {
            for (auto i = idx.first; i < idx.first + idx.second; i++) {
                if (sequences[ref_minimizers[i].seq_id].first == "S1_112901")
                    Logger::instance().debug("%lld - %d\n", it->first, ref_minimizers[i].offset);
                hits_frw.push_back((hit) {.seq_id = ref_minimizers[i].seq_id, .offset = ref_minimizers[i].offset, .hash_value = it->first});
            }
        }
        prv_hash = it->first;
    }
    sort(hits_frw.begin(), hits_frw.end());

    query_frw = hashes;
    hashes.clear();

    Logger::instance().debug("------\n");

    vector<hit> hits_rev;
    for (auto it = query_rev.begin(); it != query_rev.end(); it++) {
        auto idx = find_hit(it->first);

        if (idx.second == 0 || idx.second >= freq_th) {
            continue;
        }

        hashes.push_back(*it);

        if (it->first != prv_hash) {
            for (auto i = idx.first; i < idx.first + idx.second; i++) {
                if (sequences[ref_minimizers[i].seq_id].first == "S1_112901")
                    Logger::instance().debug("%lld - %d\n", it->first, ref_minimizers[i].offset);
                hits_rev.push_back((hit) {.seq_id = ref_minimizers[i].seq_id, .offset = ref_minimizers[i].offset, .hash_value = it->first});
            }
        }
        prv_hash = it->first;
    }
    sort(hits_rev.begin(), hits_rev.end());
    query_rev = hashes;
    hashes.clear();

    return {hits_frw, hits_rev};
}

void Sketch::get_genome_hits(string& ref_l, string& ref_r, vector<pair<uint64_t, int> >& l_frw,
                     vector<pair<uint64_t, int> >& l_rev, vector<pair<uint64_t, int> >& r_frw, vector<pair<uint64_t, int> >& r_rev,
                     pair<vector<hit>, vector<hit> >& l_hits, pair<vector<hit>, vector<hit> >& r_hits) {
    transform(ref_l.begin(), ref_l.end(), ref_l.begin(), ::toupper);
    get_query_minimizers(&ref_l[0], 0, ref_l.size(), l_frw);
    sort(l_frw.begin(), l_frw.end());

    string q = reverse_complement(ref_l);
    get_query_minimizers(&q[0], 0, q.size(), l_rev);
    sort(l_rev.begin(), l_rev.end());

    transform(ref_r.begin(), ref_r.end(), ref_r.begin(), ::toupper);
    get_query_minimizers(&ref_r[0], 0, ref_r.size(), r_frw);
    sort(r_frw.begin(), r_frw.end());

    q = reverse_complement(ref_r);
    get_query_minimizers(&q[0], 0, q.size(), r_rev);
    sort(r_rev.begin(), r_rev.end());

    Logger::instance().debug("Genomic Hits:\n");
    l_hits = get_hits(l_frw, l_rev);
    r_hits = get_hits(r_frw, r_rev);

    Logger::instance().debug("Frw Left:\n");
#ifdef DEBUG
    for (auto i = l_frw.begin(); i != l_frw.end(); i++)
        Logger::instance().debug("%lld - %d:\n", i->first, i->second);
#endif
    Logger::instance().debug("-----\n");
}

GenomeAnchorAbs Sketch::get_genome_anchor(pair<vector<hit>, vector<hit> > ref_l_hits, pair<vector<hit>, vector<hit> > ref_r_hits,
                                          id_t id, orientation_en orientation, int l_cnt, int r_cnt) {
    GenomeAnchorAbs anchor;

    if (orientation == FRW) {
        auto l_frw = equal_range(ref_l_hits.first.begin(), ref_l_hits.first.end(), id, hit());
        auto r_frw = equal_range(ref_r_hits.first.begin(), ref_r_hits.first.end(), id, hit());
        anchor.left_cnt = l_frw.second - l_frw.first;
        anchor.right_cnt = r_frw.second - r_frw.first;
        anchor.left_genome_cnt = l_cnt;
        anchor.right_genome_cnt = r_cnt;
    }
    else {
        auto l_rev = equal_range(ref_r_hits.second.begin(), ref_r_hits.second.end(), id, hit());
        auto r_rev = equal_range(ref_l_hits.second.begin(), ref_l_hits.second.end(), id, hit());
        anchor.left_cnt = l_rev.second - l_rev.first;
        anchor.right_cnt = r_rev.second - r_rev.first;
        anchor.left_genome_cnt = r_cnt;
        anchor.right_genome_cnt = l_cnt;
    }

    return anchor;
}

vector<cut> Sketch::find_cuts(bool classify, unordered_set<hash_t> &frw_minimizers, unordered_set<hash_t> &rev_minimizers) {
    vector<hit> rev_hits;
    vector<hit> frw_hits;
    rev_hits.reserve(1024);
    frw_hits.reserve(1024);

    vector<cut> cuts;

    //Find all forward hits
    for (auto it = frw_minimizers.begin(); it != frw_minimizers.end(); it++) {
        auto idx = find_hit(*it);

        if (idx.second == 0 || idx.second >= freq_th)
            continue;

        for (auto i = idx.first; i < idx.first + idx.second; i++) {
            frw_hits.push_back((hit) {.seq_id = ref_minimizers[i].seq_id, .offset = ref_minimizers[i].offset, .hash_value = *it});
        }
    }
    sort(frw_hits.begin(), frw_hits.end());

    for (auto it = rev_minimizers.begin(); it != rev_minimizers.end(); it++) {
        auto idx = find_hit(*it);

        if (idx.second == 0 || idx.second >= freq_th)
            continue;

        for (auto i = idx.first; i < idx.first + idx.second; i++) {
            rev_hits.push_back((hit) {.seq_id = ref_minimizers[i].seq_id, .offset = ref_minimizers[i].offset, .hash_value = *it});
        }
    }
    sort(rev_hits.begin(), rev_hits.end());

    //Find cuts on each picked long read
    int MIN_HITS = 0.25 * frw_minimizers.size();
    vector<cut> frw_cuts;
    for (auto it = 0; it < frw_hits.size(); it++) {
        id_t curr_id = frw_hits[it].seq_id;

        auto curr_idx = it;

        auto r = lower_bound(rev_hits.begin(), rev_hits.end(), curr_id, hit());

        if (r != rev_hits.end() && r->seq_id == curr_id)
            continue;

        auto upper = upper_bound(frw_hits.begin(), frw_hits.end(), curr_id, hit());

        hash_size_t size = upper - frw_hits.begin() - it;
        if (upper == frw_hits.end())
            size = frw_hits.size() - it;

        it += size;

        if (size < MIN_HITS)
            continue;

        cut ans = find_range(frw_hits, curr_idx, size);

        if (ans.range.start == 0 && ans.range.end == 0)
            continue;
        ans.orientation = FRW;
        ans.size = size;
        ans.seq_id = curr_id;
        //CHANGED
        if (ans.range.end >= sequences[curr_id].second)
            ans.range.end = sequences[curr_id].second - 1;

        frw_cuts.push_back(ans);
    }
//	if (frw_cuts.size() > 200)
//                return cuts;

    MIN_HITS = 0.25 * rev_minimizers.size();
    vector<cut> rev_cuts;
    for (auto it = 0; it < rev_hits.size(); it++) {
        id_t curr_id = rev_hits[it].seq_id;

        auto curr_idx = it;

        auto f = lower_bound(frw_hits.begin(), frw_hits.end(), curr_id, hit());
        if (f != frw_hits.end() && f->seq_id == curr_id)
            continue;

        auto upper = upper_bound(rev_hits.begin(), rev_hits.end(), curr_id, hit());
        hash_size_t size = upper - rev_hits.begin() - it;

        if (upper == rev_hits.end())
            size = rev_hits.size() - it;

        it += size;

        if (size < MIN_HITS)
            continue;

        cut ans = find_range(rev_hits, curr_idx, size);
        if (ans.range.start == 0 && ans.range.end == 0)
            continue;
        ans.orientation = REV;
        ans.size = size;
        ans.seq_id = curr_id;
        //CHANGED
        if (ans.range.end >= sequences[curr_id].second)
            ans.range.end = sequences[curr_id].second - 1;

        rev_cuts.push_back(ans);
    }
//	if (frw_cuts.size() + rev_cuts.size() > 200)
//                return cuts;

    if (frw_cuts.empty() && rev_cuts.empty())
        cuts =  vector<cut>();
    else if (classify) {
        classify_reads(frw_hits, frw_cuts);
        classify_reads(rev_hits, rev_cuts);

        for (int i = 0; i < frw_cuts.size(); i++) {
            cuts.push_back(frw_cuts[i]);
        }
        for (int i = 0; i < rev_cuts.size(); i++) {
            cuts.push_back(rev_cuts[i]);
        }
    }
    else {
        for (int i = 0; i < frw_cuts.size(); i++) {
            cut curr = frw_cuts[i];
            cuts.push_back(curr);
        }
    }
    return cuts;
}

float Sketch::estimate_insertion(vector<hit> hits_l, vector<hit> hits_r, vector<pair<uint64_t, int> > l_minimizers,
                                 vector<pair<uint64_t, int> > r_minimizers, id_t id, GenomeAnchorAbs anchor) {
    if (anchor.type() != FULLY_ANCHORED)
        return -1;

    auto l = equal_range(hits_l.begin(), hits_l.end(), id, hit());
    auto r = equal_range(hits_r.begin(), hits_r.end(), id, hit());

    vector<int> read_right_locs, ref_right_locs, read_left_locs, ref_left_locs;

    int cnt = 0;
    hash_t prv_hash = INT64_MAX;

    for (auto it = l.first; it != l.second - 1 && cnt < 2; it++) {
        if (it->hash_value != (it+1)->hash_value && it->hash_value != prv_hash) {
            auto left = lower_bound(l_minimizers.begin(), l_minimizers.end(), make_pair(it->hash_value, 0));
            auto right = lower_bound(r_minimizers.begin(), r_minimizers.end(), make_pair(it->hash_value, 0));
            if (right == r_minimizers.end() || right->first != it->hash_value) {
                read_left_locs.push_back(it->offset);
                ref_left_locs.push_back(left->second);
                cnt++;
            }
        }
        prv_hash = it->hash_value;
    }

    cnt = 0;
    for (auto it = r.first; it != r.second - 1 && cnt < 2; it++) {
        if (it->hash_value != (it+1)->hash_value && it->hash_value != prv_hash) {
            auto right = lower_bound(r_minimizers.begin(), r_minimizers.end(), make_pair(it->hash_value, 0));
            auto left = lower_bound(l_minimizers.begin(), l_minimizers.end(), make_pair(it->hash_value, 0));
            if (left == r_minimizers.end() || left->first != it->hash_value) {
                read_right_locs.push_back(it->offset);
                ref_right_locs.push_back(right->second);
                cnt++;
            }
        }
        prv_hash = it->hash_value;
    }

    if (read_left_locs.empty() || read_right_locs.empty()) {
        Logger::instance().debug("<No unique minimizers found>\n");
        return -1;
    }

    int ins_sum = 0;
    for (int i = 0; i < read_left_locs.size(); i++) {
        for (int j = 0; j < read_right_locs.size(); j++) {
            int tmp = abs(abs(read_right_locs[j] - read_left_locs[i]) - (genome_anchor_distance + ref_right_locs[j] - ref_left_locs[i]));
            ins_sum += tmp;
        }
    }

    int ins_avg = ins_sum/(read_left_locs.size() * read_right_locs.size());

    return ins_avg;
}

//bool genome_anchor_validation() {
//    GenomeAnchor anchor = get_genome_anchor(ref_l_hits, ref_r_hits, curr_id);
//
//    if (anchor.left_frw < REF_MIN_FRW_L && anchor.right_frw < REF_MIN_FRW_R) {
//             anchor.right_frw << "-" << anchor.left_rev << "-" << anchor.right_rev << endl;
//        continue;
//    }
//}

vector<cut> Sketch::find_cuts(bool classify, unordered_set<hash_t> &frw_minimizers, unordered_set<hash_t> &rev_minimizers, string ref_l, string ref_r) {
    vector<hit> rev_hits;
    vector<hit> frw_hits;
    rev_hits.reserve(1024);
    frw_hits.reserve(1024);

    vector<cut> cuts;

    pair<vector<hit>, vector<hit> > hits = get_hits(frw_minimizers, rev_minimizers);
    frw_hits = hits.first;
    rev_hits = hits.second;


    int MIN_HITS_FRW = 0.25 * frw_minimizers.size();
    int MIN_HITS_REV = 0.25 * rev_minimizers.size();

    vector<pair<uint64_t, int> > l_frw, l_rev, r_frw, r_rev;
    pair<vector<hit>, vector<hit> > ref_l_hits, ref_r_hits;
    get_genome_hits(ref_l, ref_r, l_frw, l_rev, r_frw, r_rev, ref_l_hits, ref_r_hits);

    Logger::instance().debug("Left frw: %d\tLeft rev: %d\tright frw: %d\tright rev: %d\n", l_frw.size(), l_rev.size(), r_frw.size(), r_rev.size());

    int REF_MIN_FRW_L = 0.25 * l_frw.size();
    int REF_MIN_FRW_R = 0.25 * r_frw.size();

    int REF_MIN_REV_L = 0.25 * l_rev.size();
    int REF_MIN_REV_R = 0.25 * r_rev.size();

    bool bimodal = false;
    bool single_peak = false;

    vector<pair<int, id_t> > estimations;

    //Find cuts on each picked long read
    float est = -1;
    float est_sum = 0;
    int est_cnt = 0;

    vector<cut> frw_cuts;
    for (auto it = 0; it < frw_hits.size(); it++) {
        id_t curr_id = frw_hits[it].seq_id;
        auto curr_idx = it;

        auto upper = upper_bound(frw_hits.begin(), frw_hits.end(), curr_id, hit());
        hash_size_t size = upper - frw_hits.begin() - curr_idx;
        if (upper == frw_hits.end())
            size = frw_hits.size() - curr_idx;

        auto rev = equal_range(rev_hits.begin(), rev_hits.end(), curr_id, hit());
        hash_size_t r_size = rev.second - rev.first;

        it += size;

        if (size < MIN_HITS_FRW || r_size >= MIN_HITS_REV) {
            continue;
        }

        Logger::instance().debug("* %s\n", sequences[curr_id].first.c_str());

        GenomeAnchorAbs anchor = get_genome_anchor(ref_l_hits, ref_r_hits, curr_id, FRW, l_frw.size(), r_frw.size());
        Logger::instance().debug("Min left anchor: %d\tMin right anchor: %d\n", anchor.left_genome_cnt, anchor.right_genome_cnt);
        Logger::instance().debug("left anchor: %d\tright anchor: %d\n", anchor.left_cnt, anchor.right_cnt);

        if (!anchor.is_anchored()) {
            continue;
        }

        cut ans = find_range(frw_hits, curr_idx, size);

        if (ans.type == BIMODAL)
            bimodal = true;
        else {
            if (anchor.type() == LEFT_ANCHORED)
                ans.type = PARTIAL_LEFT;
            else if (anchor.type() == RIGHT_ANCHORED)
                ans.type = PARTIAL_RIGHT;
        }

        float new_est = estimate_insertion(ref_l_hits.first, ref_r_hits.first, l_frw, r_frw, curr_id, anchor);
        estimations.push_back({new_est, curr_id});

        Logger::instance().debug("Insertion Estimate  = %f\n", new_est);

        if (new_est > 10) {
            est_cnt++;
            est_sum += new_est;
        }
        else if (new_est == -1) {
            if (ans.type != BIMODAL) {
                if (anchor.type() == LEFT_ANCHORED)
                    ans.type = PARTIAL_LEFT;
                else if (anchor.type() == RIGHT_ANCHORED)
                    ans.type = PARTIAL_RIGHT;
            }
        }
        else {
            Logger::instance().debug("Dropped\n", new_est);
            continue;
        }

        if (ans.range.start == 0 && ans.range.end == 0)
            continue;

        ans.orientation = FRW;
        ans.size = size;
        ans.seq_id = curr_id;

        if (ans.range.end >= sequences[curr_id].second)
            ans.range.end = sequences[curr_id].second - 1;

        ans.estimated_insertion = est;

        frw_cuts.push_back(ans);
    }

//	if (frw_cuts.size() > 200)
//                return cuts;

//========================= REV ===========================

    vector<cut> rev_cuts;
    for (auto it = 0; it < rev_hits.size(); it++) {
        id_t curr_id = rev_hits[it].seq_id;

        auto curr_idx = it;

        auto upper = upper_bound(rev_hits.begin(), rev_hits.end(), curr_id, hit());
        hash_size_t size = upper - rev_hits.begin() - it;

        if (upper == rev_hits.end())
            size = rev_hits.size() - it;

        auto frw = equal_range(frw_hits.begin(), frw_hits.end(), curr_id, hit());
        hash_size_t f_size = frw.second - frw.first;

        it += size;

        if (size < MIN_HITS_REV || f_size >= MIN_HITS_FRW) {
            continue;
        }

        Logger::instance().debug("* %s\n", sequences[curr_id].first.c_str());

        GenomeAnchorAbs anchor = get_genome_anchor(ref_l_hits, ref_r_hits, curr_id, REV, l_rev.size(), r_rev.size());
        Logger::instance().debug("Min left anchor: %d\tMin right anchor: %d\n", anchor.left_genome_cnt, anchor.right_genome_cnt);
        Logger::instance().debug("left anchor: %d\tright anchor: %d\n", anchor.left_cnt, anchor.right_cnt);

        if (!anchor.is_anchored()) {
            continue;
        }

        Logger::instance().debug("* %s\n", sequences[curr_id].first.c_str());

        cut ans = find_range(rev_hits, curr_idx, size);

        if (ans.type == BIMODAL)
            bimodal = true;

        float new_est = estimate_insertion(ref_r_hits.second, ref_l_hits.second, r_rev, l_rev, curr_id, anchor);
        estimations.push_back({new_est, curr_id});

        Logger::instance().debug("Insertion Estimate  = %f\n", new_est);

        if (new_est > 10) {
            est_cnt++;
            est_sum += new_est;
        }
        else if (new_est == -1) {
            if (ans.type != BIMODAL) {
                if (anchor.type() == LEFT_ANCHORED)
                    ans.type = PARTIAL_RIGHT;
                else if (anchor.type() == RIGHT_ANCHORED)
                    ans.type = PARTIAL_LEFT;
            }
        }
        else {
            Logger::instance().debug("Dropped\n", new_est);
            continue;
        }

        if (ans.range.start == 0 && ans.range.end == 0)
            continue;

        ans.orientation = REV;
        ans.size = size;
        ans.seq_id = curr_id;

        if (ans.range.end >= sequences[curr_id].second)
            ans.range.end = sequences[curr_id].second - 1;

        ans.estimated_insertion = est;

        rev_cuts.push_back(ans);
    }
//	if (frw_cuts.size() + rev_cuts.size() > 200)
//                return cuts;

    if (frw_cuts.empty() && rev_cuts.empty())
        cuts =  vector<cut>();
    else {
        for (int i = 0; i < frw_cuts.size(); i++) {
            cuts.push_back(frw_cuts[i]);
        }
        for (int i = 0; i < rev_cuts.size(); i++) {
            cuts.push_back(rev_cuts[i]);
        }
    }

    float est_avg;
    if (est_cnt == 0)
        est_avg = -1;
    else
        est_avg = (float)est_sum/(float)est_cnt;
    Logger::instance().debug("Average Insertion Length            : %f\n", est_avg);

    return cuts;
}