#ifndef __PROCESS_PARTITION__
#define __PROCESS_PARTITION__

#include <string>

#include "sketch.h"
#include "cut_ranges.h"
#include "progressbar.h"
#include "p3_partition.h"
#include "insertion_assembler.h"

struct cluster {
    string chrName;
    int cluster_id;
    int pt_start;
    int pt_end;
    pair<vector<read_cut_info >, classified_cuts> reads;
};

class ProcessPartition {
    private:
        int total;
        int max_len;
        FILE *fo_vcf;
        Sketch* sketch;
        FILE *fo_vcf_lq;
        int long_no = 0;
        int single_no = 0;
        int bimodal_no = 0;
        int LENFLAG	= 1000;
        int max_threads = 1;
        string reference_name;
        ProgressBar* progress;
        int processed_cnt = 0;
        cut_ranges* extractor;
        InsertionAssembler* ia;
        p3_partition* partition;
        int MAX_REF_LEN	= 300000000;

    public:
        ProcessPartition(int maxThreads, const string &lrPath, const string &datPath,
                         const string &partition, const string &range, int max_len,
                         const string &reference, const string &prefix, const string &vcf_name);
        cluster get_cluster();
        void thread_process(int tid);
        void process();
};


#endif
