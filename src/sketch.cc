#include <cmath>
#include <deque>
#include <mutex>
#include <thread>
#include <string>
#include <iostream>
#include <algorithm>
#include "logger.h"
#include "common.h"
#include "MurmurHash3.h"

#include "sketch.h"

#include <cstring>

#include <chrono>
#include <ratio>
#include <ctime>

using namespace std;

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
    kmer_size = k;
    dat_path = dp;
    window_size = w;
    load();
    compute_freq_th();
    claspChain = ClaspChain();
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

pair<vector<cut>, int> Sketch::query(vector<string> &reads, bool long_insertion, string& ref_l, string& ref_r) {
    pair<vector<cut>, int> cut_results = find_cuts_all(ref_l, ref_r, long_insertion, reads);
    return cut_results;
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

//TODO: Only use unique ones on the genome
vector<hit> Sketch::get_hits_new(vector<pair<uint64_t, int> >& query_frw) {
    vector<hit> hits_frw;
    vector<pair<hash_t, int> > hashes;

    hash_t prv_hash = UINT64_MAX;
    for (auto it = query_frw.begin(); it != query_frw.end(); it++) {
        auto idx = find_hit(it->first);

        if (idx.second == 0 || idx.second >= freq_th) {
//            if (idx.second == 0) {
            continue;
        }

        hashes.push_back(*it);

//        if (it->first != prv_hash) {
            for (auto i = idx.first; i < idx.first + idx.second; i++) {
                hits_frw.push_back((hit) {.seq_id = ref_minimizers[i].seq_id, .hash_value = it->first,
                                          .offset = ref_minimizers[i].offset, .genome_offset = (offset_t)it->second});
            }
//        }
        prv_hash = it->first;
    }
    sort(hits_frw.begin(), hits_frw.end());

    query_frw = hashes;
    hashes.clear();

    return hits_frw;
}

void Sketch::get_genome_hits_new(string& ref, vector<minimizer_t>& minimizers, vector<hit>& candidates) {
    transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
    get_query_minimizers(&ref[0], 0, ref.size(), minimizers);

    sort(minimizers.begin(), minimizers.end());

    candidates = get_hits_new(minimizers);
}

vector<seed> Sketch::create_seeds_new(vector<hit>& hits, int start, int size) {
    vector<seed> seeds;

//    ofstream fout;
//    fout.open("fragments.blast", std::ios_base::app);

    for (auto i = start; i < start + size; i++) {
        seed curr_seed;
        curr_seed.qPos = hits[i].offset;
        curr_seed.sPos = hits[i].genome_offset;
        curr_seed.len = kmer_size;
        seeds.push_back(curr_seed);
//        string ort;
//        string frg = sequences[hits[i].seq_id].first + ort + "\t" + "chrt" + "\t" + "100.0" + "\t" + to_string(kmer_size) + "\t" + "0" + "\t" + "0" +
//                     "\t" + to_string(hits[i].offset) + "\t" + to_string(hits[i].offset + kmer_size - 1) + "\t" + to_string(hits[i].genome_offset) + "\t"
//                     + to_string(hits[i].genome_offset + kmer_size - 1)
//                     + "\t" + "15.0" + "\t" + "32.0" + "\n";
//        fout << frg;
    }

//    Logger::instance().debug("SD: %d | ", seeds.size());

//    fout.close();

    return seeds;
}

MaxChainInfo Sketch::get_genome_anchor_new_new(vector<hit>& left_anchor_hits, int start, int size) {
//    ofstream fout;
//    fout.open("fragments.blast", std::ios_base::app);

//    for (auto i = left.first; i != left.second; i++) {
//        vector<int> g_offsets;
//        for (int j = 0; j < left_anchor_minimizers.size(); j++) {
//            if (left_anchor_minimizers[j].first == i->hash_value) {
//                string frg = sequences[id].first + "-frw-l\t" + "chrt" + "\t" + "100.0" + "\t" + "15" + "\t" + "0" + "\t" + "0" +
//                             "\t" + to_string(i->offset) + "\t" + to_string(i->offset + 14) + "\t" + to_string(left_anchor_minimizers[j].second) + "\t"
//                             + to_string(left_anchor_minimizers[j].second + 14)
//                             + "\t" + "15.0" + "\t" + "32.0" + "\n";
//                fout << frg;
//            }
//        }
//    }

    auto t1 = chrono::high_resolution_clock::now();
    vector<seed> lseeds = create_seeds_new(left_anchor_hits, start, size);
    auto t21 = chrono::high_resolution_clock::now();
    seed_time += std::chrono::duration<double, std::milli>(t21-t1).count();
    auto t11 = chrono::high_resolution_clock::now();
    MaxChainInfo max_chain_l = claspChain.get_max_chain(lseeds);
    auto t2 = chrono::high_resolution_clock::now();
    chaining_time += std::chrono::duration<double, std::milli>(t2-t1).count();
    clasp_time += std::chrono::duration<double, std::milli>(t2-t11).count();
    int number_of_minimizers = (max_chain_l.gaps_size + max_chain_l.score)/kmer_size;
    Logger::instance().debug("C: (Q: %d-%d, G: %d-%d, L: %d, S: %f, K: %d) | ", max_chain_l.qrange.first,
                             max_chain_l.qrange.second, max_chain_l.rrange.first, max_chain_l.rrange.second, max_chain_l.len,
                             max_chain_l.score, number_of_minimizers);

//    for (auto i = right.first; i != right.second; i++) {
//        vector<int> g_offsets;
//        for (int j = 0; j < right_anchor_minimizers.size(); j++) {
//            if (right_anchor_minimizers[j].first == i->hash_value) {
//                string frg = sequences[id].first + "-frw-r\t" + "chrt" + "\t" + "100.0" + "\t" + "15" + "\t" + "0" + "\t" + "0" +
//                             "\t" + to_string(i->offset) + "\t" + to_string(i->offset + 14) + "\t" + to_string(right_anchor_minimizers[j].second) + "\t"
//                             + to_string(right_anchor_minimizers[j].second + 14)
//                             + "\t" + "15.0" + "\t" + "32.0" + "\n";
//                fout << frg;
//            }
//        }
//    }

//    fout.close();

    return max_chain_l;
}

inline int Sketch::max_kmer_count(MaxChainInfo chain, int anchor_length, int read_length, int genome_minimizer_count) {
    int qdistance = chain.qrange.second - chain.qrange.first;
    int gdistance = chain.rrange.second - chain.rrange.first;
    float read_distance_coefficient = (float)qdistance/(float)gdistance;
    int p_left = 1;
    if (chain.qrange.first < chain.rrange.first * read_distance_coefficient) {
        p_left = (int)((chain.rrange.first - chain.qrange.first) / read_distance_coefficient);
    }

    int p_right = anchor_length;
    if (read_length - chain.qrange.second < (anchor_length - chain.rrange.second) * read_distance_coefficient) {
        p_right = chain.rrange.second + (int)((read_length - chain.qrange.second)/read_distance_coefficient);
    }

    Logger::instance().debug("PL: %d | PR: %d | RL: %d | ", p_left, p_right, read_length);

    return ((p_right - p_left)/anchor_length) * genome_minimizer_count;
}

void Sketch::find_left_cuts(vector<hit>& read_candidates, vector<cut>& candidates_cut_info, int anchor_length, orientation_en orientation,
                            int genome_minimizers_cnt, unordered_set<id_t>& insertion_candidates) {
    Logger::instance().debug("&&Minimizers: %d\n", genome_minimizers_cnt);
    vector<int> left_reads_remaining_size, right_reads_remaining_size;

    for (auto it = 0; it < read_candidates.size(); it++) {
        id_t curr_id = read_candidates[it].seq_id;
        auto curr_idx = it;

        auto upper = upper_bound(read_candidates.begin(), read_candidates.end(), curr_id, hit());
        hash_size_t size = upper - read_candidates.begin() - curr_idx;
        if (upper == read_candidates.end())
            size = read_candidates.size() - curr_idx;

        it += size;

        if (insertion_candidates.size() != 0) {
            if (insertion_candidates.find(curr_id) == insertion_candidates.end()) {
//                Logger::instance().debug("%-30s DROPPED(IM)\n", sequences[curr_id].first.c_str());
                continue;
            }
        }

        int minimizer_cutoff = (int)((1.0 * GENOME_ANCHOR_CUTOFF/anchor_length) * genome_minimizers_cnt);
//        if (size < max(3, minimizer_cutoff)) {
        if (size < 3) {
//            tmp_cnt++;
            continue;
        }

        Logger::instance().debug("%-30s", sequences[curr_id].first.c_str());
        //move curr_id to first argument and drop orientation
        MaxChainInfo max_chain = get_genome_anchor_new_new(read_candidates, curr_idx, size);
        Logger::instance().debug(" S: %d | ", size);
        tmp_cnt++;

//        S2_228529                     C: (Q: 16-116, G: 487-587, L: 101, S: 90.000000, K: 6) | DROPPED(LM)

        int number_of_minimizers = (max_chain.gaps_size + max_chain.score)/kmer_size;
        if (number_of_minimizers < 0.25 * max_kmer_count(max_chain, anchor_length, (int)sequences[curr_id].second, genome_minimizers_cnt)) {
            Logger::instance().debug("DROPPED(LM)\n");
            continue;
        }

//        if (anchor_length - max_chain.rrange.second > BOUNDARY_DISTANCE_CUTOFF) {
//            Logger::instance().debug("DROPPED(BR)\n");
//            continue;
//        }
        if (max_chain.qrange.first > BOUNDARY_DISTANCE_CUTOFF && max_chain.rrange.first > BOUNDARY_DISTANCE_CUTOFF) {
            Logger::instance().debug("DROPPED(QS)\n");
            continue;
        }
        if (max_chain.qrange.second - max_chain.qrange.first < GENOME_ANCHOR_CUTOFF) {
            Logger::instance().debug("DROPPED(QL)\n");
            continue;
        }
        int qdistance = max_chain.qrange.second - max_chain.qrange.first;
        int gdistance = max_chain.rrange.second - max_chain.rrange.first;
        float read_distance_coefficient = (float)qdistance/(float)gdistance;
        if (qdistance - gdistance > 0.15 * gdistance) {
            Logger::instance().debug("DROPPED(TH)\n");
            continue;
        }

        int pivot = max(anchor_length - GENOME_ANCHOR_LEN, max_chain.rrange.first);
        int skip = (pivot - max_chain.rrange.second) * read_distance_coefficient;

        if (max_chain.qrange.second + skip > (int)sequences[curr_id].second - 10) {
            Logger::instance().debug("DROPPED(NER)\n");
            continue;
        }

        cut read_cut_info;
        read_cut_info.seq_id = curr_id;
        read_cut_info.orientation = orientation;
        read_cut_info.range.start = max(0, max_chain.qrange.second + skip);
        read_cut_info.range.end = sequences[curr_id].second - 1;
        read_cut_info.genome_range.start = max_chain.rrange.first;
        read_cut_info.genome_range.end = pivot;
        read_cut_info.breakpoint_distance = abs(pivot - anchor_length);

        Logger::instance().debug("Q: [%4d - %4d, %4d] | G: [%4d - %4d] | P: %d | ++\n", read_cut_info.range.start, read_cut_info.range.end,
                                 read_cut_info.range.end - read_cut_info.range.start, read_cut_info.genome_range.start,
                                 read_cut_info.genome_range.end, read_cut_info.breakpoint_distance);

        candidates_cut_info.push_back(read_cut_info);
    }
}

void Sketch::find_right_cuts(vector<hit>& read_candidates, vector<cut>& candidates_cut_info, int anchor_length, orientation_en orientation,
                             int genome_minimizers_cnt, unordered_set<id_t>& insertion_candidates) {
    Logger::instance().debug("&&Minimizers: %d\n", genome_minimizers_cnt);
    vector<int> left_reads_remaining_size, right_reads_remaining_size;
    int sum = 0;
    int cnt = 0;

    for (auto it = 0; it < read_candidates.size(); it++) {
        int begin = 0, end = 0;
        id_t curr_id = read_candidates[it].seq_id;
        auto curr_idx = it;

        auto upper = upper_bound(read_candidates.begin(), read_candidates.end(), curr_id, hit());
        hash_size_t size = upper - read_candidates.begin() - curr_idx;
        if (upper == read_candidates.end())
            size = read_candidates.size() - curr_idx;

        it += size;

        if (insertion_candidates.size() != 0) {
            if (insertion_candidates.find(curr_id) == insertion_candidates.end()) {
//                Logger::instance().debug("%-30s DROPPED(IM)\n", sequences[curr_id].first.c_str());
                continue;
            }
        }

        int minimizer_cutoff = (int)((1.0 * GENOME_ANCHOR_CUTOFF/anchor_length) * genome_minimizers_cnt);
//        if (size < max(3, minimizer_cutoff)) {
        if (size < 3) {
//            tmp_cnt++;
            continue;
        }

        Logger::instance().debug("%-30s", sequences[curr_id].first.c_str());
        //move curr_id to first argument and drop orientation
        MaxChainInfo max_chain = get_genome_anchor_new_new(read_candidates, curr_idx, size);
        Logger::instance().debug(" S: %d | ", size);
        int qdistance = max_chain.qrange.second - max_chain.qrange.first;
        int gdistance = max_chain.rrange.second - max_chain.rrange.first;
        tmp_cnt++;

        //TODO: add to paper
        int number_of_minimizers = (max_chain.gaps_size + max_chain.score)/kmer_size;
        if (number_of_minimizers < 0.25 * max_kmer_count(max_chain, anchor_length, (int)sequences[curr_id].second, genome_minimizers_cnt)) {
            Logger::instance().debug("DROPPED(LM)\n");
            continue;
        }

//        if (max_chain.rrange.first > BOUNDARY_DISTANCE_CUTOFF) {
//            Logger::instance().debug("DROPPED(BR)\n ");
//            continue;
//        }
        if (sequences[curr_id].second - max_chain.qrange.second > BOUNDARY_DISTANCE_CUTOFF &&
            anchor_length - max_chain.rrange.second > BOUNDARY_DISTANCE_CUTOFF) {
            Logger::instance().debug("DROPPED(QE)\n");
            continue;
        }
        if (max_chain.qrange.second - max_chain.qrange.first < GENOME_ANCHOR_CUTOFF) {
            Logger::instance().debug("DROPPED(QL)\n");
            continue;
        }

        if (qdistance - gdistance > 0.15 * gdistance) {
            Logger::instance().debug("DROPPED(TH)\n");
            continue;
        }

        int pivot = min(GENOME_ANCHOR_LEN, max_chain.rrange.second);
        int skip = (pivot - max_chain.rrange.first) * (float)qdistance/(float)gdistance;

        if (max_chain.qrange.first + skip < 10) {
            Logger::instance().debug("DROPPED(NEL)\n");
            continue;
        }

        cut read_cut_info;
        read_cut_info.seq_id = curr_id;
        read_cut_info.orientation = orientation;
        read_cut_info.range.end = min(sequences[curr_id].second - 1, max(max_chain.qrange.first + skip, 0));
        read_cut_info.range.start = 0;
        read_cut_info.genome_range.start = pivot;
        read_cut_info.genome_range.end = max_chain.rrange.second;
        read_cut_info.breakpoint_distance = pivot;

        Logger::instance().debug("Q: [%4d - %4d, %4d] | G: [%4d - %4d] | P: %d | **\n", read_cut_info.range.start, read_cut_info.range.end,
                                 read_cut_info.range.end - read_cut_info.range.start, read_cut_info.genome_range.start,
                                 read_cut_info.genome_range.end, read_cut_info.breakpoint_distance);

        candidates_cut_info.push_back(read_cut_info);
    }
}

void Sketch::merge_candidates(vector<cut>& left_candidates, vector<cut>& right_candidates, vector<cut>& merged_candidates,
                              cut_stats& stats, bool long_insertion) {
    int i = 0, j = 0;
    cut dummy;
    dummy.seq_id = UINT32_MAX;
    left_candidates.push_back(dummy);
    right_candidates.push_back(dummy);

    while (i < left_candidates.size() || j < right_candidates.size()) {
         if (left_candidates[i].seq_id < right_candidates[j].seq_id) {
            left_candidates[i].type = PARTIAL_LEFT;
            merged_candidates.push_back(left_candidates[i]);
            stats.left_cuts_size.push_back(left_candidates[i].range.end - left_candidates[i].range.start);
            i++;
        }
        else if (right_candidates[j].seq_id < left_candidates[i].seq_id) {
            right_candidates[j].type = PARTIAL_RIGHT;
            merged_candidates.push_back(right_candidates[j]);
            stats.right_cuts_size.push_back(right_candidates[j].range.end - right_candidates[j].range.start);
            j++;
        }
        else {
            if (left_candidates[i].seq_id == UINT32_MAX)
                break;
            if (left_candidates[i].range.start > right_candidates[j].range.end) {
//                int read_start_distance = left_candidates[i].range.start - right_candidates[j].range.end;
//                int genome_left_remaining = abs(left_candidates[i].genome_range.start - right_candidates[j].genome_range.end);
//                int genome_right_remaining = abs(right_candidates[j].genome_range.start - left_candidates[i].genome_range.end);
//                if (read_start_distance <= genome_left_remaining && read_start_distance <= genome_right_remaining) {
                    cut merged_cut = left_candidates[i];
                    merged_cut.range.start = right_candidates[j].range.end;
                    merged_cut.range.end = left_candidates[i].range.start;

                    if (merged_cut.range.end - merged_cut.range.start < left_candidates[i].breakpoint_distance +
                        right_candidates[j].breakpoint_distance) {
                        int extension = ((left_candidates[i].breakpoint_distance + right_candidates[j].breakpoint_distance)
                                - (merged_cut.range.end - merged_cut.range.start)) / 2;
                        merged_cut.range.start = max(0, merged_cut.range.start - extension);
                        merged_cut.range.end = min((int)sequences[merged_cut.seq_id].second, merged_cut.range.end + extension);
                    }

                    merged_cut.type = OVERLAPPING_READ;
                    Logger::instance().debug("Setting Overlap\n");
                    merged_cut.estimated_insertion = -1;
                    merged_candidates.push_back(merged_cut);
                    stats.overlapping_cnt++;
//                }
//                else {
//                    Logger::instance().debug("%-30s LS > RS | No decision\n", sequences[left_candidates[i].seq_id].first.c_str());
//                }
                i++;
                j++;
            }
            else {
                cut merged_cut = left_candidates[i];
                merged_cut.range.end = right_candidates[j].range.end;
                merged_cut.estimated_insertion = abs(merged_cut.range.end - merged_cut.range.start -
                        (left_candidates[i].breakpoint_distance + right_candidates[j].breakpoint_distance));
//                merged_cut.estimated_insertion = abs(merged_cut.range.end - merged_cut.range.start -
//                                                     (left_candidates[i].genome_range.start - right_candidates[j].genome_range.end));
                if (merged_cut.range.end - merged_cut.range.start < left_candidates[i].breakpoint_distance +
                                                                    right_candidates[j].breakpoint_distance) {
                    int extension = ((left_candidates[i].breakpoint_distance + right_candidates[j].breakpoint_distance)
                                     - (merged_cut.range.end - merged_cut.range.start)) / 2;
                    //TODO: Fix extensions
//                    int left_extension = left_candidates[i].breakpoint_distance - (merged_cut.range.end - merged_cut.range.start);
//                    int right_extension = right_candidates[j].breakpoint_distance - (merged_cut.range.end - merged_cut.range.start);
                    merged_cut.range.start = max(0, merged_cut.range.start - extension);
                    merged_cut.range.end = min((int)sequences[merged_cut.seq_id].second, merged_cut.range.end + extension);
                }
                merged_cut.type = BIMODAL;
                i++;
                j++;
                if (!long_insertion && merged_cut.estimated_insertion < 10) {
                    Logger::instance().debug("DROPPED(%s, I=%d, R=%d-%d)\n", sequences[merged_cut.seq_id].first.c_str(),
                                             merged_cut.estimated_insertion, merged_cut.range.start, merged_cut.range.end);
                    continue;
                }
                stats.bimodal_sum += merged_cut.estimated_insertion;
                stats.bimodal_cnt++;
                merged_candidates.push_back(merged_cut);
            }
        }
    }
    sort(stats.left_cuts_size.begin(), stats.left_cuts_size.end());
    sort(stats.right_cuts_size.begin(), stats.right_cuts_size.end());
    Logger::instance().debug("%d\n", merged_candidates.size());
}

void Sketch::get_unique_minimizers(vector<string> &reads, unordered_set<hash_t>& insertion_minimizers) {
    vector<minimizer_t> minimizers_vec;
    for (int i = 0; i < reads.size(); i++) {
        transform(reads[i].begin(), reads[i].end(), reads[i].begin(), ::toupper);
        get_query_minimizers(&reads[i][0], i, reads[i].size(), minimizers_vec);
    }
    for (int i = 0; i < minimizers_vec.size(); i++) {
        insertion_minimizers.insert(minimizers_vec[i].first);
//        Logger::instance().debug("%lld\n", minimizers_vec[i].first);
    }
}

void Sketch::get_insertion_minimizers_new(vector<string>& short_reads, string& genome_left, string& genome_right, unordered_set<id_t>& candidates) {
    unordered_set<hash_t> reads_minimizers;
//    Logger::instance().debug("Reads: \n");
    get_unique_minimizers(short_reads, reads_minimizers);
//    for (auto it = reads_minimizers.begin(); it != reads_minimizers.end(); it++) {
//        Logger::instance().debug("%llu\n", *it);
//    }

    unordered_set<hash_t> genome_minimizers;

    transform(genome_left.begin(), genome_left.end(), genome_left.begin(), ::toupper);
    vector<string> genome_segments;
    genome_segments.push_back(genome_left);
    genome_segments.push_back(genome_right);
//    Logger::instance().debug("Genome: \n");
    get_unique_minimizers(genome_segments, genome_minimizers);
//    for (auto it = genome_minimizers.begin(); it != genome_minimizers.end(); it++) {
//        Logger::instance().debug("%llu\n", *it);
//    }

    vector<hit> candidates_vec;

    int insertion_minimizers_size = 0;

    Logger::instance().debug("Insertion Minimizers:\n");
    for (auto it = reads_minimizers.begin(); it != reads_minimizers.end(); it++) {
        if (genome_minimizers.find(*it) == genome_minimizers.end()) {
//            Logger::instance().debug("%llu\n", *it);
            insertion_minimizers_size++;
            auto idx = find_hit(*it);
            for (auto i = idx.first; i < idx.first + idx.second; i++) {
                hit tmp_hit;
                tmp_hit.seq_id = ref_minimizers[i].seq_id;
                tmp_hit.hash_value = *it;
                candidates_vec.push_back(tmp_hit);
            }
        }
    }

    Logger::instance().debug("Insertion Minimizer Size: %d\n", insertion_minimizers_size);

    sort(candidates_vec.begin(), candidates_vec.end());

    id_t prv_id = UINT32_MAX;
    hash_t prv_hash = UINT64_MAX;
    int cnt = 0;
    int cutoff = INSERTION_MINIMIZERS_CUTOFF * insertion_minimizers_size;
    for (int i = 0; i < candidates_vec.size(); i++) {
        if (candidates_vec[i].seq_id == prv_id) {
            if (candidates_vec[i].hash_value != candidates_vec[i-1].hash_value) {
                cnt++;
            }
        }
        else {
            if (cnt >= cutoff) {
                candidates.insert(prv_id);
            }
            if (prv_id != UINT32_MAX)
                Logger::instance().debug("%s(%d) ", sequences[prv_id].first.c_str(), cnt);
            prv_id = candidates_vec[i].seq_id;
            cnt = 1;
        }
    }
    Logger::instance().debug("\n");
}

vector<cut> Sketch::find_cuts_with_chain(string ref_l, string ref_r, cut_stats& stats, orientation_en orientation,
                                         unordered_set<id_t>& insertion_candidates, bool long_insertion) {
    vector<minimizer_t> left_frw_minimizers, left_rc_minimizers, right_frw_minimizers, right_rc_minimizers;
    vector<hit> left_frw_candidates, left_rc_candidates, right_frw_candidates, right_rc_candidates;

    auto t1 = chrono::high_resolution_clock::now();
    get_genome_hits_new(ref_l, left_frw_minimizers, left_frw_candidates);
    get_genome_hits_new(ref_r, right_frw_minimizers, right_frw_candidates);
    auto t2 = chrono::high_resolution_clock::now();
    hits_time += std::chrono::duration<double, std::milli>(t2-t1).count();

    Logger::instance().debug("Minimizers Size: %d\n", left_frw_minimizers.size());
    Logger::instance().debug("Hits Size: %d\n", left_frw_candidates.size());

    vector<cut> left_frw_cuts, right_frw_cuts;

    t1 = chrono::high_resolution_clock::now();
    Logger::instance().debug("--- LEFT ---\n");
    find_left_cuts(left_frw_candidates, left_frw_cuts, ref_l.size(), orientation, left_frw_minimizers.size(), insertion_candidates);
    Logger::instance().debug("--- RIGHT ---\n");
    find_right_cuts(right_frw_candidates, right_frw_cuts, ref_r.size(), orientation, right_frw_minimizers.size(), insertion_candidates);
    t2 = chrono::high_resolution_clock::now();
    finding_time += std::chrono::duration<double, std::milli>(t2-t1).count();

    vector<cut> final_cuts;
    t1 = chrono::high_resolution_clock::now();
    merge_candidates(left_frw_cuts, right_frw_cuts, final_cuts, stats, long_insertion);
    t2 = chrono::high_resolution_clock::now();
    merging_time += std::chrono::duration<double, std::milli>(t2-t1).count();
    Logger::instance().debug("%d\n", final_cuts.size());

    #ifdef DEBUG
    sort(final_cuts.begin(), final_cuts.end());
    for (int  i = 0; i < final_cuts.size(); i++) {
        Logger::instance().debug("%-30s: %4d-%4d (%4d) | ", sequences[final_cuts[i].seq_id].first.c_str(), final_cuts[i].range.start,
                                 final_cuts[i].range.end, final_cuts[i].range.end - final_cuts[i].range.start);
        Logger::instance().debug("I: %d | ", final_cuts[i].estimated_insertion);

        if (final_cuts[i].type == BIMODAL)
            Logger::instance().debug("T: B\n");
        else if (final_cuts[i].type == PARTIAL_LEFT)
            Logger::instance().debug("T: L\n");
        else if (final_cuts[i].type == PARTIAL_RIGHT)
            Logger::instance().debug("T: R\n");
        else if (final_cuts[i].type == OVERLAPPING_READ)
            Logger::instance().debug("T: O\n");
        else
            Logger::instance().debug("T: M\n");
    }
#endif

    return final_cuts;
}

void inline readjust(cut& read_cut_info, int left_extension, int right_extension, bool reorient) {
    if (read_cut_info.type == PARTIAL_LEFT && left_extension != -1) {
        read_cut_info.range.end = min((int)read_cut_info.range.end, read_cut_info.range.start + left_extension);
    }
    else if (read_cut_info.type == PARTIAL_RIGHT && right_extension != -1) {
        read_cut_info.range.start = max((int)read_cut_info.range.start, read_cut_info.range.end - right_extension);
    }
    if (reorient) {
        if (read_cut_info.type == PARTIAL_LEFT)
            read_cut_info.type = PARTIAL_RIGHT;
        else if (read_cut_info.type == PARTIAL_RIGHT)
            read_cut_info.type = PARTIAL_LEFT;
    }
}

pair<vector<cut>, int> Sketch::find_cuts_all(string& ref_l, string& ref_r, bool long_insertion, vector<string>& short_reads) {
    vector<cut> cuts, cuts_2;

    unordered_set<id_t> insertion_candidates_frw, insertion_candidates_rc;
    if (!short_reads.empty()) {
        int read_len = short_reads[0].size();
        int left_cut_size = min(read_len, (int)ref_l.size());
        int right_cut_size = min(read_len, (int)ref_r.size());

        string genome_left = ref_l.substr(ref_l.size() - left_cut_size, left_cut_size);
        string genome_left_rc = reverse_complement(genome_left);

        string genome_right = ref_r.substr(0, right_cut_size);
        string genome_right_rc = reverse_complement(genome_right);

        vector<string> short_reads_rc(short_reads.size());
        for (int i = 0; i < short_reads.size(); i++) {
            short_reads_rc.push_back(reverse_complement(short_reads[i]));
        }

        get_insertion_minimizers_new(short_reads, genome_left, genome_right, insertion_candidates_frw);
        get_insertion_minimizers_new(short_reads_rc, genome_left_rc, genome_right_rc, insertion_candidates_rc);

//        Logger::instance().debug("Insertion Minimizers: F(%d), R(%d)\n", insertion_candidates_frw.size(), insertion_candidates_rc.size());
    }

    auto t1 = chrono::high_resolution_clock::now();
    cut_stats frw_stats, rc_stats;
    Logger::instance().debug("--- FRW ---\n");
    cuts = find_cuts_with_chain(ref_l, ref_r, frw_stats, FRW, insertion_candidates_frw, long_insertion);
    Logger::instance().debug("--- REV ---\n");
    cuts_2 = find_cuts_with_chain(reverse_complement(ref_r), reverse_complement(ref_l), rc_stats, REV, insertion_candidates_rc, long_insertion);
    auto t2 = chrono::high_resolution_clock::now();

    int insertion_estimation = -1;
    t1 = chrono::high_resolution_clock::now();
    int left_extension = -1, right_extension = -1;
    if (frw_stats.bimodal_cnt + rc_stats.bimodal_cnt > 2) {
        insertion_estimation = (frw_stats.bimodal_sum + rc_stats.bimodal_sum) / (frw_stats.bimodal_cnt + rc_stats.bimodal_cnt);
        left_extension = insertion_estimation + GENOME_ANCHOR_LEN;
        right_extension = insertion_estimation + GENOME_ANCHOR_LEN;
        Logger::instance().debug("IE: %d | LE: %d | RE: %d\n", insertion_estimation, left_extension, right_extension);
    }
    else {
        vector<int> left_sizes(frw_stats.left_cuts_size.size() + rc_stats.right_cuts_size.size());
        vector<int> right_sizes(frw_stats.right_cuts_size.size() + rc_stats.left_cuts_size.size());

        std::merge(frw_stats.left_cuts_size.begin(), frw_stats.left_cuts_size.end(), rc_stats.right_cuts_size.begin(),
              rc_stats.right_cuts_size.end(), left_sizes.begin());
        Logger::instance().debug("Left: \n");
        for (int i = 0; i < left_sizes.size(); i++) {
            Logger::instance().debug("%d, ", left_sizes[i]);
        }
        Logger::instance().debug("LF: \n");
        for (int i = 0; i < frw_stats.left_cuts_size.size(); i++) {
            Logger::instance().debug("%d, ", frw_stats.left_cuts_size[i]);
        }
        Logger::instance().debug("RR: \n");
        for (int i = 0; i < rc_stats.right_cuts_size.size(); i++) {
            Logger::instance().debug("%d, ", rc_stats.right_cuts_size[i]);
        }
        Logger::instance().debug("\n");
        std::merge(frw_stats.right_cuts_size.begin(), frw_stats.right_cuts_size.end(), rc_stats.left_cuts_size.begin(),
              rc_stats.left_cuts_size.end(), right_sizes.begin());
        Logger::instance().debug("Right: \n");
        for (int i = 0; i < right_sizes.size(); i++) {
            Logger::instance().debug("%d, ", right_sizes[i]);
        }
        Logger::instance().debug("\n");

        if (left_sizes.size() > 7) {
            left_extension = left_sizes[left_sizes.size() - 4];
        }
        else if (left_sizes.size() > 2) {
            left_extension = left_sizes[left_sizes.size()/2];
        }
        if (right_sizes.size() > 7) {
            right_extension = right_sizes[right_sizes.size() - 4];
        }
        else if (right_sizes.size() > 2) {
            right_extension = right_sizes[right_sizes.size()/2];
        }
        Logger::instance().debug("LE: %d | RE: %d\n", left_extension, right_extension);
    }
    for (int i = 0; i < cuts.size(); i++) {
        readjust(cuts[i], left_extension, right_extension, false);
    }
    for (int i = 0; i < cuts_2.size(); i++) {
        readjust(cuts_2[i], right_extension, left_extension, true);
    }
    t2 = chrono::high_resolution_clock::now();

    cuts.insert(cuts.end(), cuts_2.begin(), cuts_2.end());

    t1 = chrono::high_resolution_clock::now();
#ifdef DEBUG
    sort(cuts.begin(), cuts.end());
    for (int  i = 0; i < cuts.size(); i++) {
        Logger::instance().debug("%-30s: %4d-%4d (%4d) | ", sequences[cuts[i].seq_id].first.c_str(), cuts[i].range.start,
                                 cuts[i].range.end, cuts[i].range.end - cuts[i].range.start);
        Logger::instance().debug("I: %d | ", cuts[i].estimated_insertion);

        if (cuts[i].type == BIMODAL)
            Logger::instance().debug("T: B\n");
        else if (cuts[i].type == PARTIAL_LEFT)
            Logger::instance().debug("T: L\n");
        else if (cuts[i].type == PARTIAL_RIGHT)
            Logger::instance().debug("T: R\n");
        else if (cuts[i].type == OVERLAPPING_READ)
            Logger::instance().debug("T: O\n");
        else
            Logger::instance().debug("T: M\n");
    }
#endif
    t2 = chrono::high_resolution_clock::now();
//    cerr << "printing: " << std::chrono::duration<double, std::milli>(t2-t1).count() << endl;

    if (frw_stats.overlapping_cnt + rc_stats.overlapping_cnt > 2) {
        return {cuts, -1};
    }
    if (insertion_estimation < 300 && insertion_estimation > 0) {
        return {cuts, -(frw_stats.bimodal_cnt + rc_stats.bimodal_cnt)};
    }
    return {cuts, frw_stats.bimodal_cnt + rc_stats.bimodal_cnt};
}