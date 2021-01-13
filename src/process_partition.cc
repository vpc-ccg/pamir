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

    sketch = new Sketch(datPath);
    extractor = new cut_ranges(lrPath, true);
    ia = new InsertionAssembler(sketch, extractor);

    progress = new ProgressBar(80);
    char comment[20];
    sprintf(comment, "%10d / %-10d", 0, total);
    progress->update((0.0/(float)total) * 100, comment);

    string out_vcf = prefix + "/" + vcf_name + ".vcf";
    fo_vcf 				= fopen(out_vcf.c_str(), "w");
    string out_vcf_lq = prefix  + "/" + vcf_name + "_LOW_QUAL.vcf";
    fo_vcf_lq 			= fopen(out_vcf_lq.c_str(), "w");
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
    if (ans.reads.first.size() != 0) {
        processed_cnt += 1;
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
    char tmp[2000];

    spoa::Graph graph{};
    auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 2, -32, -64, -1);
    genome reference(reference_name.c_str());
    aligner al(max_len + 2010);

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
        logger_out += " + Reference             : " + ref_part + "\n";

        // if the genomic region is too big
        if (ref_end - ref_start > MAX_REF_LEN)
            continue;

        // holding the calls info, can be used to detect the repeated calls, etc.
        vector<tuple<string, int, int, string, int, float> > reports;
        vector<tuple<string, int, int, string, int, float> > reports_lq;

        vector<pair<string, int> > consensus;

        //BIMODAL
        if (p.reads.second.bimodal) {
            logger_out += " + Type                  : bimodal(left: " + to_string(p.reads.second.left_cuts.size()) +
                    ", right: " + to_string(p.reads.second.right_cuts.size()) +
                    ", bimodal: " + to_string(p.reads.second.bimodal_cuts.size()) + ")\n\n";

            graph.Clear();
            for (int i = 0; i < p.reads.second.bimodal_cuts.size(); i++) {
                auto alignment = alignment_engine->Align(p.reads.second.bimodal_cuts[i].first.second, graph);
                graph.AddAlignment(alignment, p.reads.second.bimodal_cuts[i].first.second);
            }
            for (int i = 0; i < p.reads.second.left_cuts.size(); i++) {
                auto alignment = alignment_engine->Align(p.reads.second.left_cuts[i].first.second, graph);
                graph.AddAlignment(alignment, p.reads.second.left_cuts[i].first.second);
            }
            for (int i = 0; i < p.reads.second.right_cuts.size(); i++) {
                auto alignment = alignment_engine->Align(p.reads.second.right_cuts[i].first.second, graph);
                graph.AddAlignment(alignment, p.reads.second.right_cuts[i].first.second);
            }

            auto msa = graph.GenerateMultipleSequenceAlignment(true);

            pair<string, pair<int, int>> ans = cut_consensus_bimodal(msa, p.reads.second.left_cuts.size(), p.reads.second.right_cuts.size(),
                                                                     p.reads.second.bimodal_cuts.size());

            consensus.push_back({ans.first, p.reads.second.left_cuts.size() + p.reads.second.bimodal_cuts.size() + p.reads.second.right_cuts.size()});

            bimodal_no++;
        }
        else {
            //SINGLE PEAK
            if (p.reads.second.right_cuts.size() < 2) {
                logger_out += " + Type              : single peak\n\n";

                graph.Clear();
                for (int i = 0; i < p.reads.second.left_cuts.size(); i++) {
                    auto alignment = alignment_engine->Align(p.reads.second.left_cuts[i].first.second, graph);
                    graph.AddAlignment(alignment, p.reads.second.left_cuts[i].first.second);
                }

                auto msa = graph.GenerateMultipleSequenceAlignment(true);

                pair<string, pair<int, int>> ans = cut_consensus_single(msa);

                consensus.push_back({ans.first, p.reads.second.left_cuts.size()});

                single_no++;
            }
            //LONG INSERTION
            else if (!p.reads.second.left_cuts.empty() && !p.reads.second.right_cuts.empty()) {
                logger_out += " + Type              : long insertion\n\n";

                vector<string> reads_1, reads_2;
                vector<string> left_reads, right_reads;

                graph.Clear();
                for (int i = 0; i < p.reads.second.left_cuts.size(); i++) {
                    auto alignment = alignment_engine->Align(p.reads.second.left_cuts[i].first.second, graph);
                    graph.AddAlignment(alignment, p.reads.second.left_cuts[i].first.second);
                    reads_1.push_back(p.reads.second.left_cuts[i].first.second);
                }
                auto msa_1 = graph.GenerateMultipleSequenceAlignment(true);
                string cons_1 = msa_1[msa_1.size() - 1];
                cons_1.erase(std::remove(cons_1.begin(), cons_1.end(), '-'), cons_1.end());

                graph.Clear();
                for (int i = 0; i < p.reads.second.right_cuts.size(); i++) {
                    auto alignment = alignment_engine->Align(p.reads.second.right_cuts[i].first.second, graph);
                    graph.AddAlignment(alignment, p.reads.second.right_cuts[i].first.second);
                    reads_2.push_back(p.reads.second.right_cuts[i].first.second);
                }
                auto msa_2 = graph.GenerateMultipleSequenceAlignment(true);
                string cons_2 = msa_2[msa_2.size() - 1];
                cons_2.erase(std::remove(cons_2.begin(), cons_2.end(), '-'), cons_2.end());

                string left, right;
                al.align(ref_part, cons_1);

                //left flank
                if (al.get_left_anchor() > 5) {
                    left = cons_1;
                    right = cons_2;
                    left_reads = reads_1;
                    right_reads = reads_2;
                }
                    //right flank
                else if (al.get_right_anchor() > 5) {
                    left = cons_2;
                    right = cons_1;
                    left_reads = reads_2;
                    right_reads = reads_1;
                }

                pair<string, int> ans = ia->assemble(left_reads, right_reads);

                consensus.push_back(ans);

                long_no++;
            }
        }

        for (int i = 0; i < consensus.size(); i++) {
            int contig_support = consensus[i].second;

            logger_out += "\n\n>>>>> Length: " + to_string(consensus[i].first.size()) +
                    " Support: " + to_string(contig_support) + " Contig: " + consensus[i].first + "\n";

            //Align consensus
            al.align(ref_part, consensus[i].first);
            string ex_msg = "";
            if (al.extract_calls(cluster_id, reports_lq, reports, contig_support, ref_start, ">>>", ex_msg) == 0) {
                string rc_contig = reverse_complement(consensus[i].first);
                al.align(ref_part, rc_contig);
                al.extract_calls(cluster_id, reports_lq, reports, contig_support, ref_start, "<<<", ex_msg);
            }

            logger_out += ex_msg;
        }

        if (log_buffer % 1000 == 0) {
            log_mutex.lock();
            Logger::instance().info(logger_out.c_str());
            log_mutex.unlock();
            logger_out.clear();
            log_buffer = 0;
        }

        //print
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

        append_vcf(chrName, out_ref, reports, p.cluster_id, vcf_info, vcf_info_del);
        n_buffer++;
        if (0 == n_buffer % 500) {
            vcf_mutex.lock();
            fprintf(fo_vcf, "%s", vcf_info.c_str());
            n_buffer = 0;
            vcf_info.clear();
            vcf_mutex.unlock();
        }
        append_vcf(chrName, out_ref_lq, reports_lq, p.cluster_id, vcf_info_lq, vcf_info_del);
        n_buffer2++;
        if (n_buffer2 == 0)
            n_buffer2++;
        if (0 == n_buffer2 % 500) {
            vcf_mutex.lock();
            fprintf(fo_vcf_lq, "%s", vcf_info_lq.c_str());
            n_buffer2 = 0;
            vcf_info_lq.clear();
            vcf_mutex.unlock();
        }
    }

    if (logger_out.size() > 0) {
        log_mutex.lock();
        Logger::instance().info(logger_out.c_str());
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
}
