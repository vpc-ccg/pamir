#include "process_partition.h"

#include <tuple>
#include <mutex>
#include <thread>
#include <iterator>
#include <algorithm>

#include "common.h"
#include "genome.h"
#include "aligner.h"
#include "cut_ranges.h"
#include "spoa/spoa.hpp"

using namespace std;

mutex pt_mutex, log_mutex, vcf_mutex;

ProcessPartition::ProcessPartition(int maxThreads, const string &lrPath, const string &datPath,
                                   const string &partition_file, const string &range, int max_len,
                                   const string &reference_name, const string &prefix, const string &vcf_name) :
                                   max_threads(maxThreads), max_len(max_len), reference_name(reference_name) {
    partition = new p3_partition(partition_file, range);
    total = partition->get_total();

    progress = new ProgressBar(80);
    char comment[20];
    sprintf(comment, "%10d / %-10d", 0, total);
    progress->update((0.0/(float)total) * 100, comment);

    sketch = new Sketch(datPath, 0, 0);
    extractor = new cut_ranges(lrPath, true);
    ia = new InsertionAssembler(sketch, extractor);

    string out_vcf = prefix + "/" + vcf_name + ".vcf";
    fo_vcf 				= fopen(out_vcf.c_str(), "w");
    string out_vcf_lq = prefix  + "/" + vcf_name + "_LOW_QUAL.vcf";
    fo_vcf_lq 			= fopen(out_vcf_lq.c_str(), "w");

    string log_file = prefix + "/" + "all" + ".log";
    fo_log = fopen(log_file.c_str(), "w");

    fprintf(fo_vcf, "%s", "##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description=\"All filters passed\">\n##INFO=<ID=Cluster,Number=1,Type=Integer,Description=\"ID of the cluster the variant is extracted from\">\n##INFO=<ID=Support,Number=1,Type=Integer,Description=\"Number of reads/contigs supporting the contig\">\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n");
}

inline string cluster_header(int cluster_id, int sr_size) {
    string header = "";
    header += "-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n";
    header += " + Cluster ID            : " + to_string(cluster_id) + "\n";
    header += " + Short Reads Count     : " + to_string(sr_size) + "\n";
    return header;
}

cluster ProcessPartition::get_cluster() {
    cluster ans;
    char comment[20];
    pt_mutex.lock();
    ans.reads = partition->read_partition();
    ans.chrName = partition->get_reference();
    ans.pt_end = partition->get_end();
    ans.pt_start = partition->get_start();
    ans.cluster_id = partition->get_id();
    ans.estimated_insertion = partition->get_estimated_insertion();
    ans.cluster_type = ans.reads.second.cluster_type;
    if (ans.reads.first.size() != 0) {
        processed_cnt += 1;
        cerr << "cluster: " << processed_cnt << endl;
        Logger::instance().debug("+ Cluster ID            : %d\n", processed_cnt);
        sprintf(comment, "%10d / %-10d", processed_cnt, total);
        progress->update(((float)processed_cnt/(float)total) * 100, comment);
    }
    pt_mutex.unlock();
    return ans;
}

void ProcessPartition::thread_process(int tid) {
    string out_ref; out_ref.reserve(4);
    string out_ref_lq; out_ref_lq.reserve(4);
    string vcf_info = ""; vcf_info.reserve(10000000);
    string vcf_info_lq = ""; vcf_info_lq.reserve(10000000);
    string vcf_info_del = ""; vcf_info.reserve(10000000);
    string logger_out = ""; logger_out.reserve(10000000);
    int n_buffer         =   0;
    int n_buffer2         =   0;
    int log_buffer = 0;

    spoa::Graph graph{};
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 10, -2, -15, -7);
//    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 2, -32, -64, -1);
    genome reference(reference_name.c_str());
    aligner al(max_len + 2010);

    string cluster_type = "";

    while (1) {
        cluster p = get_cluster();

        // end of the partition file
        if (!p.reads.first.size())
            break;

        log_buffer++;

        // cluster has too many or too few reads
        if ( p.reads.first.size() > 7000 || p.reads.first.size() <= 2 ) {
            logger_out += cluster_header(p.cluster_id, p.reads.first.size());
            logger_out += "INFO: Skipped Processing - Too few or Too many short reads\n";
            continue;
        }
        else if (p.reads.second.size == 0) {
            logger_out += cluster_header(p.cluster_id, p.reads.first.size());
            logger_out += "INFO: Skipped Processing - No long reads found\n";
            continue;
        }
        else if ( p.reads.second.size > 150 || p.reads.second.size <= 2 ) {
            logger_out += cluster_header(p.cluster_id, p.reads.first.size());
            logger_out += " + Long Reads Count      : " + to_string(p.reads.second.size) + "\n";
            logger_out += "INFO: Skipped Processing - Too few or Too many long reads\n";
            continue;
        }

        string chrName  = p.chrName;
        int cluster_id  = p.cluster_id;
        int pt_start    = p.pt_start;
        int pt_end      = p.pt_end;
        int ref_start   = pt_start - LENFLAG;
        int ref_end     = pt_end   + LENFLAG;
        string ref_part = reference.get_bases_at(chrName, ref_start, ref_end);
        logger_out += cluster_header(p.cluster_id, p.reads.first.size());
        logger_out += " + Long Reads Count      : " + to_string(p.reads.second.size) + "\n";
        logger_out += " + Spanning Range        : " + chrName + ":" + to_string(pt_start) + "-" + to_string(pt_end) + "\n";
        logger_out += " + Discovery Range       : " + chrName + ":" + to_string(ref_start) + "-" + to_string(ref_end) + "\n";
        logger_out += " + Estimated Insertion   : " + to_string(p.estimated_insertion) + "\n";
        logger_out += " + Reference             : " + ref_part + "\n";

        // if the genomic region is too big
        if (ref_end - ref_start > MAX_REF_LEN)
            continue;

        // holding the calls info, can be used to detect the repeated calls, etc.
        vector<tuple<string, int, int, string, int, float> > reports;
        vector<tuple<string, int, int, string, int, float> > reports_lq;

        vector<pair<string, int> > consensus;

        //BIMODAL
        if (p.cluster_type == BIMODAL && p.reads.second.bimodal_cuts.size() > 2) {
            cluster_type = "bimodal";
            logger_out += " + Type                  : bimodal(left: " + to_string(p.reads.second.left_cuts.size()) +
                    ", right: " + to_string(p.reads.second.right_cuts.size()) +
                    ", bimodal: " + to_string(p.reads.second.bimodal_cuts.size()) + ")\n\n";

            int b_size = 0, l_size = 0, r_size = 0;

            graph.Clear();

            int l = 0, r = 0;
            for (int i = 0; i < min(7, (int)p.reads.second.bimodal_cuts.size()); i++) {
                auto alignment = alignment_engine->Align(p.reads.second.bimodal_cuts[i].sequence, graph);
                graph.AddAlignment(alignment, p.reads.second.bimodal_cuts[i].sequence);
            }
            if (p.estimated_insertion > 0) {
                for (int i = 0; i < min(7, (int)p.reads.second.left_cuts.size()); i++) {
                    auto alignment = alignment_engine->Align(p.reads.second.left_cuts[i].sequence, graph);
                    graph.AddAlignment(alignment, p.reads.second.left_cuts[i].sequence);
                }
                for (int i = 0; i < min(7, (int)p.reads.second.right_cuts.size()); i++) {
                    auto alignment = alignment_engine->Align(p.reads.second.right_cuts[i].sequence, graph);
                    graph.AddAlignment(alignment, p.reads.second.right_cuts[i].sequence);
                }
                l = min(7, (int)p.reads.second.left_cuts.size());
                r = min(7, (int)p.reads.second.left_cuts.size());
            }

            auto msa = graph.GenerateMultipleSequenceAlignment(true);

            pair<string, pair<int, int>> ans = cut_consensus_bimodal(msa, l, r,
                                                                     min(7, (int)p.reads.second.bimodal_cuts.size()));

            consensus.push_back({ans.first, p.reads.second.left_cuts.size() + p.reads.second.bimodal_cuts.size() + p.reads.second.right_cuts.size()});

            #ifdef DEBUG
            logger_out += "\nMSA:\n";
            for (int i = 0; i < msa.size(); i++) {
                string tmp = msa[i].insert(ans.second.first, "|");
                string tmp2 = tmp.insert(ans.second.second + 1, "|");
                logger_out += tmp2 + "\n";
            }
            #endif

            bimodal_no++;
        }
        //LONG INSERTION
        else if (!p.reads.second.left_cuts.empty() && !p.reads.second.right_cuts.empty()) {
            cluster_type = "long";
            logger_out += " + Type              : long insertion\n\n";

            vector<string> left_reads, right_reads;

            for (int i = 0; i < min((int)p.reads.second.left_cuts.size(), 7); i++) {
                left_reads.push_back(p.reads.second.left_cuts[i].sequence);
            }
            for (int i = 0; i < min((int)p.reads.second.right_cuts.size(), 7); i++) {
                right_reads.push_back(p.reads.second.right_cuts[i].sequence);
            }

            pair<string, int> ans = ia->assemble(left_reads, right_reads, logger_out);

            consensus.push_back(ans);

            long_no++;
        }

        for (int i = 0; i < consensus.size(); i++) {
            int contig_support = consensus[i].second;

            logger_out += "\n\n>>>>> Length: " + to_string(consensus[i].first.size()) +
                    " Support: " + to_string(contig_support) + " Contig: " + consensus[i].first + "\n";

            if (consensus[i].first.size() > 30000) {
                logger_out += "\nContig too long, skipping.\n";
                continue;
            }

            //Align consensus
            al.align(ref_part, consensus[i].first);
            string ex_msg = "";
            al.extract_calls(cluster_id, reports_lq, reports, contig_support, ref_start, ">>>", ex_msg);
            logger_out += ex_msg;
        }

        if (log_buffer % 1000 == 0) {
            log_mutex.lock();
            fprintf(fo_log, "%s", logger_out.c_str());
            log_mutex.unlock();
            logger_out.clear();
            log_buffer = 0;
        }

        out_ref.clear();
        for (int j = 0; j < reports.size(); j++) {
            int tmp_end = get<1>(reports[j]);
            out_ref += reference.get_base_at(chrName, tmp_end);
        }

        out_ref_lq.clear();
        for (int j =0; j < reports_lq.size(); j++) {
            int tmp_end = get<1>(reports_lq[j]);
            out_ref_lq += reference.get_base_at(chrName, tmp_end);
        }

        append_vcf_hybrid(chrName, out_ref, reports, p.cluster_id, cluster_type, p.reads.second.bimodal_cuts.size(),
                          p.reads.second.left_cuts.size(), p.reads.second.right_cuts.size(),
                          p.reads.second.misc_cuts.size(),
                          p.estimated_insertion, vcf_info, vcf_info_del);

        n_buffer++;
        if (0 == n_buffer % 500) {
            vcf_mutex.lock();
            fprintf(fo_vcf, "%s", vcf_info.c_str());
            vcf_mutex.unlock();
            n_buffer = 0;
            vcf_info.clear();
        }
        append_vcf_hybrid(chrName, out_ref_lq, reports_lq, p.cluster_id, cluster_type,
                          p.reads.second.bimodal_cuts.size(),
                          p.reads.second.left_cuts.size(), p.reads.second.right_cuts.size(),
                          p.reads.second.misc_cuts.size(),
                          p.estimated_insertion, vcf_info_lq, vcf_info_del);

        n_buffer2++;
        if (n_buffer2 == 0)
            n_buffer2++;
        if (0 == n_buffer2 % 500) {
            vcf_mutex.lock();
            fprintf(fo_vcf_lq, "%s", vcf_info_lq.c_str());
            vcf_mutex.unlock();
            n_buffer2 = 0;
            vcf_info_lq.clear();
        }
    }

    if (logger_out.size() > 0) {
        log_mutex.lock();
        fprintf(fo_log, "%s", logger_out.c_str());
        log_mutex.unlock();
    }
    if (vcf_info.size() > 0) {
        vcf_mutex.lock();
        fprintf(fo_vcf, "%s", vcf_info.c_str());
        n_buffer = 0;
        vcf_info.clear();
        vcf_mutex.unlock();
    }
    if (vcf_info.size() > 0) {
        vcf_mutex.lock();
        fprintf(fo_vcf_lq, "%s", vcf_info_lq.c_str());
        n_buffer2 = 0;
        vcf_info_lq.clear();
        vcf_mutex.unlock();
    }
}

void ProcessPartition::process() {
    thread threads[max_threads];
    for (int i = 0; i < max_threads; i++)
        threads[i] = thread(&ProcessPartition::thread_process, this, i);
    for (int i = 0; i < max_threads; i++)
        threads[i].join();
    fclose(fo_vcf);
    fclose(fo_vcf_lq);
    fclose(fo_log);
}
