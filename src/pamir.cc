#include <algorithm>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <cstdio>
#include <string>
#include <map>
#include <tuple>
#include <utility>
#include <set>
#include <vector>
#include <cmath>
#include <zlib.h>
#include "common.h"
#include "partition.h"
#include "assembler.h"
#include "genome.h"
#include "aligner.h"
#include "extractor.h"
#include "record.h"
#include "sam_parser.h"
#include "process_range.h"
#include "process_orphans.h"
#include "recalibrate.h"
#include "smoother.h"
#include "sketch.h"
#include "p2_partition.h"
#include "p3_partition.h"
#include "cut_ranges.h"
#include "insertion_assembler.h"
#include "progressbar.h"
#include "spoa/spoa.hpp"
#include "process_partition.h"
#include "partition_hybrid.h"
#include <chrono>


using namespace std;
//TODO Implement repeaet master after genome file is converted to util

/****************************************************************/
// For outputing specific log
void log_idx (const string &log_file )
{
    FILE *fin = fopen(log_file.c_str(), "rb");
    FILE *fidx = fopen((log_file + ".idx").c_str(), "wb");
    char *readline = (char*)malloc(MAX_CHAR);
    char *token = (char*)malloc(100);
    size_t idx_pos = ftell(fin);
    int l_id, offset;
    int num_inserted = 0; // to resolve skipping partition issue

    fwrite( &idx_pos, 1, sizeof(size_t), fidx); // initialize an log for partition id ZERO
    while( NULL != fgets( readline, MAX_CHAR, fin ) )
    {
        if ( 0 == strncmp("-<=*=>-*-<", readline, 10) )
        {
            fgets( readline, MAX_CHAR, fin);
            sscanf(readline, "%s %s %s %s %d\n", token, token, token, token, &l_id);
            while( l_id > num_inserted +1)
            {
                //fprintf( stdout, "size\t%d\t%lu->%s\n", num_inserted, idx_pos, readline);
                fwrite( &idx_pos, 1, sizeof(size_t), fidx);
                num_inserted++;
            }
            //fprintf( stdout, "size\t%d\t%lu->%s\n", num_inserted, idx_pos, readline);
            fwrite( &idx_pos, 1, sizeof(size_t), fidx);
            num_inserted++;
        }
        idx_pos = ftell(fin);
    }
    fclose(fin);
    fclose(fidx);
    free(readline);
}
/****************************************************************/
// Output Log from x to y-1. To output t, specify t-t+1
int output_log (const string &log_file, const string &range)
{
    static unsigned int start = -1, end = -1;
    static vector<size_t> offsets;
    if (start == -1) {
        char *dup = strdup(range.c_str());
        char *tok = strtok(dup, "-");
        if (!tok) start = 0;
        else {
            start = atol(tok), tok = strtok(0, "-");
            end = tok ? atol(tok) : start+1;
        }
        free(dup);
        //free(tok);
        fprintf(stdout, "extraction [%u, %u]\n", start, end-1);

        FILE *fidx = fopen((log_file + ".idx").c_str(), "rb");
        size_t offset;
        while (fread(&offset, 1, sizeof(size_t), fidx) == sizeof(size_t))
            offsets.push_back(offset);
        fclose(fidx);
    }

    FILE *fi, *fo, *foidx;
    int sz, i;
    int cluster_id;
    int num_cluster = 0, num_read = 0;
    const int MAXB = 8096;
    char *pref = (char*)malloc(MAX_CHAR);
//	char name[MAXB], read[MAXB];
    string c_file = range + ".log";
    fo = fopen(c_file.c_str(), "w");
    fclose(fo);
    reset:
//	assert(start < offsets.size());
    if (start >= offsets.size() || start >= end)
        return 0;
    //fprintf(stderr,"Seeking to %d--%d (%lu)\n", start, end, offsets[start]);

    fi = fopen(log_file.c_str(), "rb");
    fo = fopen(c_file.c_str(), "a");
    fseek(fi, offsets[start++], SEEK_SET);
    fgets(pref, MAXB, fi);
    if ( 0 != strncmp("-<=*=>-*-<", pref, 10) )
    {	exit(1); fprintf(stderr, "Incorrect Start at %s", pref);
    }
    fprintf( fo, "%s", pref);

    fgets(pref, MAXB, fi);
    while ( 0 != strncmp("-<=*=>-*-<", pref, 10) )
    {
        fprintf( fo, "%s", pref);
        fgets(pref, MAXB, fi);
    }

    num_read = 0;

    fclose(fi);
    fclose(fo);
    delete pref;
    if ( num_read == 0)
        goto reset;
    return num_cluster;
}
/*******************************************************************/
void print_header(const string &header_file, const string &reference)
{
    string header_info;
    header_info.reserve(16384);

    // access genome information
    genome toread(reference.c_str());
    char *absref = new char[1000];
    char *baseref = new char[1000];
    strcpy(absref,reference.c_str());
    baseref = strtok(absref,"/");
    char *prevref = new char[500];
    while(baseref!=NULL)
    {
        strcpy(prevref,baseref);
        baseref=strtok(NULL,"/");
    }
    toread.load_next();

    header_info =	"##fileformat=VCFv4.2\n"
                     "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
                     "##reference=";
    header_info+=   prevref;
    header_info+=	"\n##source=Pamir\n"
                     "##ALT=<ID=<INS>,Type=String,Description=\"Insertion of novel sequence\">\n"
                     "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
                     "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n"
                     "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End coordinate of this variant\">\n"
                     "##INFO=<ID=Cluster,Number=1,Type=Integer,Description=\"ID of the cluster the variant is extracted from\">\n"
                     "##INFO=<ID=Support,Number=1,Type=Integer,Description=\"Number of reads/contigs supporting the contig\">\n"
                     "##INFO=<ID=SEQ,Number=1,Type=String,Description=\"Variant sequence\">\n";

    string prevName="";
    string name = toread.get_name();
    int ssize = toread.get_size();
    while (name!=prevName && ssize!=0)
    {
        header_info +=	"##contig=<ID=";	header_info+=	name;
        header_info	+=	",length=";		 	header_info+=	to_string(ssize);
        header_info	+=	">\n";
        prevName=name;
        toread.load_next();
        name  =toread.get_name();
        ssize = toread.get_size();
    }
    header_info +=	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    FILE *fo = fopen(header_file.c_str(),"w");
    fprintf(fo, "%s", header_info.c_str() );
    fclose(fo);
}
/****************************************************************/
void assemble (const string &partition_file, const string &reference, const string &range, const string &name, int max_len, int read_length, const string &prefix)
{
    const double MAX_AT_GC 		= 0.7;
    const int MAX_REF_LEN		= 300000000;
    int LENFLAG					= 1000;//500;//1000;
    char *line 					= new char[MAX_CHAR];
    string out_vcf = prefix + "/" + name + ".vcf";
    FILE *fo_vcf 				= fopen(out_vcf.c_str(), "w");
//	string out_vcf_del = out_vcf.substr(0,out_vcf.rfind("."))+string("_DELS")+out_vcf.substr(out_vcf.rfind("x")+1,out_vcf.length());
//	string out_vcf_lq = out_vcf.substr(0,out_vcf.rfind("."))+string("_LOWQUAL")+out_vcf.substr(out_vcf.rfind("x")+1,out_vcf.length());
    string out_vcf_del = prefix + "/" + name + "_DELS.vcf";
    string out_vcf_lq = prefix  + "/" + name + "_LOW_QUAL.vcf";
    FILE *fo_vcf_del 			= fopen(out_vcf_del.c_str(), "w");
    FILE *fo_vcf_lq 			= fopen(out_vcf_lq.c_str(), "w");

    assembler as(max_len, 15);
    genome ref(reference.c_str());
    map<string,string> chroms;
    genome_partition pt(partition_file, range);
    aligner al(max_len + 2010 );

    string tmp_ref; tmp_ref.reserve(4);
    string tmp_ref_lq; tmp_ref_lq.reserve(4);
    string vcf_info=""; vcf_info.reserve(10000000);
    string vcf_info_lq=""; vcf_info_lq.reserve(10000000);
    string vcf_info_del=""; vcf_info.reserve(10000000);
    const int MAX_BUFFER = 500;
    int n_buffer         =   0;
    int n_buffer2         =   0;

    while (1)
    {
        auto p 			= pt.read_partition();
        // end of the partition file
        if ( !p.size() )
            break;

        // cluster has too many or too few reads
        if ( p.size() > 7000 || p.size() <= 2 ) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID      : %d\n", pt.get_id());
            Logger::instance().info(" + Reads Count     : %lu\n", p.size());
            Logger::instance().info("INFO: Skipped Processing - Too few or Too many reads\n");
            continue;
        }
        string chrName  = pt.get_reference();
        int cluster_id  = pt.get_id();
        int pt_start    = pt.get_start();
        int pt_end      = pt.get_end();
        int ref_start   = pt_start - LENFLAG;
        int ref_end     = pt_end   + LENFLAG;
        string ref_part = ref.get_bases_at(chrName, ref_start, ref_end);
        Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
        Logger::instance().info(" + Cluster ID      : %d\n", cluster_id);
        Logger::instance().info(" + Reads Count     : %lu\n", p.size());
        Logger::instance().info(" + Spanning Range  : %s:%d-%d\n", chrName.c_str(), pt_start, pt_end);
        Logger::instance().info(" + Discovery Range : %s:%d-%d\n", chrName.c_str(), ref_start, ref_end);
        Logger::instance().info(" + Reference       : %s\n\n", ref_part.c_str());
        // if the genomic region is too big
        if (ref_end - ref_start > MAX_REF_LEN)
            continue;

        // holding the calls info, can be used to detect the repeated calls, etc.
        vector< tuple< string, int, int, string, int, float > > reports;//reports.clear();
        vector< tuple< string, int, int, string, int, float > > reports_lq;//reports.clear();

        vector<string> reads;
        for (int i =0;i<p.size();i++)
        {
            reads.push_back(p[i].first.second);
        }
        auto contigs    = as.assemble(reads);
        for ( auto &contig: contigs )
        {
            int contig_support		= contig.read_information.size();
            int con_len 			= contig.data.length();
            if( check_AT_GC(contig.data, MAX_AT_GC) == 0 || (con_len <= read_length && contig_support <= 1) || con_len > max_len + 400 ) continue;

            Logger::instance().info("\n\n>>>>> Length: %d Support: %d Contig: %s\n", con_len, contig_support, contig.data.c_str());
            for(int z=0;z<contig.read_information.size();z++)
                Logger::instance().info("%s %s %d %d\n",
                                        contig.read_information[z].name.c_str(),
                                        contig.read_information[z].seq.c_str(),
                                        contig.read_information[z].location,
                                        contig.read_information[z].location_in_contig);
            al.align(ref_part, contig.data);
            if(al.extract_calls(cluster_id, reports_lq, reports, contig_support, ref_start,">>>")==0)
            {
                string rc_contig = reverse_complement(contig.data);
                al.align(ref_part, rc_contig);
                al.extract_calls(cluster_id, reports_lq, reports, contig_support, ref_start, "<<<");
            }
        }
        //print_calls new version
        tmp_ref.clear();//string tmp_ref = "";
        for (int j =0; j <reports.size();j++)
        {
            int tmp_end = get<1>(reports[j]);
            tmp_ref += ref.get_base_at(chrName, tmp_end);
            //tmp_ref += ref.extract(chrName, tmp_end, tmp_end);
        }

        tmp_ref_lq.clear();//string tmp_ref = "";
        for (int j =0; j <reports_lq.size();j++)
        {
            int tmp_end = get<1>(reports_lq[j]);
            tmp_ref_lq += ref.get_base_at(chrName, tmp_end);
        }
        append_vcf( chrName, tmp_ref, reports, pt.get_id(), vcf_info, vcf_info_del);
        n_buffer++;
        if ( 0 == n_buffer%MAX_BUFFER )
        {
            fprintf( fo_vcf, "%s", vcf_info.c_str());
            n_buffer = 0;
            vcf_info.clear();
        }
        append_vcf( chrName, tmp_ref_lq, reports_lq, pt.get_id(), vcf_info_lq, vcf_info_del);
        n_buffer2++;
        /*if (n_buffer==0)
            n_buffer++;
        if ( 0 == n_buffer%MAX_BUFFER )
        {
            fprintf( fo_vcf_del, "%s", vcf_info_del.c_str());
            n_buffer = 0;
            vcf_info_del.clear();
        }*/
        if(n_buffer2 ==0)
            n_buffer2++;
        if ( 0 == n_buffer2%MAX_BUFFER )
        {
            fprintf( fo_vcf_lq, "%s", vcf_info_lq.c_str());
            n_buffer2 = 0;
            vcf_info_lq.clear();
        }
    }
    // Sanity check for the last record
    if ( 0 < vcf_info.size()){	fprintf( fo_vcf, "%s", vcf_info.c_str());}
    //if ( 0 < vcf_info_del.size()){	fprintf( fo_vcf_del, "%s", vcf_info_del.c_str());}
    if ( 0 < vcf_info_lq.size()){	fprintf( fo_vcf_lq, "%s", vcf_info_lq.c_str());}
    fclose(fo_vcf);
    fclose(fo_vcf_lq);
    fclose(fo_vcf_del);
}
/*********************************************************************************************/
void sketch (const string &longread, const string &dat_path, int k, int w)
{
    Sketch lr_sketch = Sketch(longread, dat_path, k, w);
}
/*********************************************************************************************/
void find_reads (const string &partition_file, const string &dat_path, const string &range, const string &ref, const int len)
{
    genome_partition pt(partition_file, range);

    ProgressBar progress(80);
    char comment[20];
    int cnt = 0;
    int total = pt.get_total();
    sprintf(comment, "%10d / %-10d", cnt, total);
    progress.update(((float)cnt/(float)total) * 100, comment);

    int genome_anchor_len = 600;
//    int genome_anchor_len = 550;
    int genome_anchor_distance = 650;
//    int genome_anchor_len = 3 * len;
//    int genome_anchor_distance = genome_anchor_len + 2 * len;

    Sketch lr_sketch = Sketch(dat_path, len, genome_anchor_distance);

    //TODO FIX NAME - Get ouput as an argument
    string p2 = "partition-p2-" + range;
    p2_partition pt_2(p2, true);

    vector<pair<id_t, range_s> > ranges;
    ranges.reserve(1000000);

    genome reference(ref.c_str());

    while (1)
    {
        cnt++;
        auto p 			= pt.read_partition();

        // end of the partition file
        if ( !p.size() )
            break;

        sprintf(comment, "%10d / %-10d", cnt, total);
        progress.update(((float)cnt/(float)total) * 100, comment);

        pair<vector<cut>, int> cut_results;
        vector<cut> cut_candidates;

        string chrName  = pt.get_reference();
        int pt_start    = pt.get_start();
        int pt_end      = pt.get_end();

        Logger::instance().debug("------------------------------------\n");
        Logger::instance().debug("+ Cluster ID            : %d\n", cnt);
        Logger::instance().debug("+ Range                 : %d-%d\n", pt_start, pt_end);
        // cluster has too many or too few reads
        if ( p.size() > 7000 || p.size() <= 2 ) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID      : %d\n", pt.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.size());
            Logger::instance().info("INFO: Skipped Processing - Too few or Too many short reads\n");

            pt_2.add_cuts(p, cut_candidates, pt_start, pt_end, chrName, -1);
            continue;
        }

        vector<string> reads;
        reads.reserve(p.size());
        for (int i = 0; i < p.size(); i++) {
            reads.push_back(p[i].first.second);
        }

        int breakpoint = pt_start +  (pt_end - pt_start)/2;
        int left_start = breakpoint - genome_anchor_len;
        int left_end = breakpoint;
        int right_start = breakpoint;
        int right_end = breakpoint + genome_anchor_len;

        string ref_l = reference.get_bases_at(chrName, left_start, left_end);
        string ref_r = reference.get_bases_at(chrName, right_start, right_end);

        cut_results = lr_sketch.query(reads, false, ref_l, ref_r);
        cut_candidates = cut_results.first;

        int estimated_insertion = cut_results.second;
        for (int i = 0; i < cut_candidates.size(); i++) {
            ranges.push_back({cut_candidates[i].seq_id, cut_candidates[i].range});
//            if (cut_candidates[i].estimated_insertion != -1 && cut_candidates[i].estimated_insertion < estimated_insertion)
//                estimated_insertion = cut_candidates[i].estimated_insertion;
        }

        if (cut_candidates.size() == 0) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID            : %d\n", pt.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.size());
            Logger::instance().info("INFO: Skipped Processing - No Long reads found\n");

            pt_2.add_cuts(p, cut_candidates, pt_start, pt_end, chrName, -1);
            continue;
        }

        pt_2.add_cuts(p, cut_candidates, pt_start, pt_end, chrName, estimated_insertion);
    }

    sort(ranges.begin(), ranges.end(), [](const pair<id_t, range_s> &a, const pair<id_t, range_s> &b) {
        return a.first < b.first;
    });
    uint32_t size = ranges.size();
    uint32_t negs = 0;

    if (!ranges.empty()) {
        for (uint32_t i = 1; i < ranges.size(); i++) {
            if (ranges[i].first != ranges[i-1].first)
                continue;
            else {
                ranges[i].second.start = min(ranges[i].second.start, ranges[i-1].second.start);
                ranges[i].second.end = max(ranges[i].second.end, ranges[i-1].second.end);
                ranges[i-1].first = -1;
                negs++;
            }
        }
    }

    sort(ranges.begin(), ranges.end(), [](const pair<id_t, range_s> &a, const pair<id_t, range_s> &b) {
        return a.first < b.first;
    });

    size -= negs;

    ofstream fout("extract_ranges.dat", ios::out | ios::binary);
    fout.write((char*)&size, sizeof(uint32_t));
    fout.write((char*)&ranges[0], size * sizeof(pair<id_t, range_s>));
    fout.close();

    cerr << "Hits Time: " << lr_sketch.hits_time << endl;
    cerr << "Find Time: " << lr_sketch.finding_time << endl;
    cerr << "Merging Time: " << lr_sketch.merging_time << endl;
    cerr << "Chaining Time: " << lr_sketch.chaining_time << endl;
    cerr << "Chain Sorting Time: " << lr_sketch.claspChain.sort_time << endl;
    cerr << "Clasp Time: " << lr_sketch.clasp_time << endl;
    cerr << "Seed Time: " << lr_sketch.seed_time << endl;
    cerr << "Building Chain: " << lr_sketch.claspChain.chain_time << endl;
    cerr << "Cnt: " << lr_sketch.tmp_cnt << endl;
}
/*********************************************************************************************/
void extract_reads(const string &partition_file, const string &longread, const string &range, const string &p3_name,
                   const string& dat_path, const string& ranges_file)
{
    p2_partition pt_2(partition_file, range);

    ProgressBar progress(80);
    char comment[20];
    int cnt = 0;
    int total = pt_2.get_total();
    sprintf(comment, "%10d / %-10d", cnt, total);
    progress.update(((float)cnt/(float)total) * 100, comment);
    cut_ranges ranges = cut_ranges(longread, dat_path, ranges_file);

    ranges.extract();

    //TODO FIX NAME + get output name as argument
    string p3 = "partition-p3-" + p3_name;
    p3_partition pt_3(p3, true);

    aligner al(30000);

    while (1) {
        cnt++;

        auto p 			= pt_2.read_partition();
        if ( !p.first.size() )
            break;

        sprintf(comment, "%10d / %-10d", cnt, total);
        progress.update(((float)cnt/(float)total) * 100, comment);

        vector<p3_read_s> cuts;

        string chrName  = pt_2.get_reference();
        int pt_start    = pt_2.get_start();
        int pt_end      = pt_2.get_end();

        // cluster has too many or too few reads
        if ( p.first.size() > 7000 || p.first.size() <= 2 ) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID      : %d\n", pt_2.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.first.size());
            Logger::instance().info("INFO: Skipped Processing - Too few or Too many short reads\n");

            pt_3.add_reads(p.first, cuts, pt_start, pt_end, chrName, pt_2.get_estimated_insertion());
            continue;
        }
            // cluster has too many or too few reads
        else if ( p.second.size() == 0 ) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID      : %d\n", pt_2.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.first.size());
            Logger::instance().info("INFO: Skipped Processing - No long reads found\n");

            pt_3.add_reads(p.first, cuts, pt_start, pt_end, chrName, pt_2.get_estimated_insertion());
            continue;
        }

        for (int i = 0; i < p.second.size(); i++) {
            pair<string, string> tmp = ranges.get_cut(p.second[i].id, p.second[i].range.start,
                                                      p.second[i].range.end);
            string cut;
            if (p.second[i].orientation == REV) {
                cut = reverse_complement(tmp.second);
            }
            else
                cut = tmp.second;
            cuts.push_back((p3_read_s) {.name = tmp.first, .sequence = cut, .range = p.second[i].range,
                    .type = p.second[i].type, .orientation = p.second[i].orientation});
        }

        pt_3.add_reads(p.first, cuts, pt_start, pt_end, chrName, pt_2.get_estimated_insertion());
    }
}
/*********************************************************************************************/
void consensus (const string &partition_file, const string &reference, const string lr_path, const string dat_path,
                const string &range, const string &name, int max_len, const string &prefix)
{
    ProcessPartition processor = ProcessPartition(8, lr_path, dat_path, partition_file, range, max_len, reference,
                                                  prefix, name);
    processor.process();
}
/*********************************************************************************************/
int main(int argc, char **argv)
{
    try {

        if (argc < 2) throw "Usage:\tpamir [mode=(?)]";

        string mode = argv[1];
        if (mode == "remove_concordant") {
            if (argc != 8) throw "Usage:\tpamir remove_concordant [sam-file] [output] [output type extension] oea? orphan? matched_ratio";
            extractor ext(argv[2], argv[3], argv[4], atoi(argv[5]), atoi(argv[6]), stod(argv[7]));
        }
        else if (mode == "partition") {
            if ( 8 == argc)
            {
                string log_path = argv[3];
                log_path+=".log";
                Logger::instance().info.set_file(log_path.c_str());
//                Logger::instance().info.set_prefix("[PAMIR] ").toggle_time();
                genome_partition pt(atoi(argv[4]), argv[2],  argv[5], argv[6], argv[7], argv[3]) ;
                pt.cluster_reads();
            }
            else{ throw "Usage:\tpamir partition [read-file] [output-file] [threshold] [ [orphan-contig] [oea2orphan] ] [mate_file]"; }
        }
        else if (mode == "get_cluster") {
            if (argc != 4) throw "Usage:\tpamir get_cluster [partition-file] [range]";
//            Logger::instance().info.set_prefix("[PAMIR] ").toggle_time();
            string log_path = string(argv[3])+".cluster";
            Logger::instance().info.set_file(log_path.c_str());
            genome_partition pt(argv[2], argv[3], true);
            pt.output_partitions();
        }
        else if (mode=="header"){
            if (argc !=4) throw "Usage:3 parameters needed\tpamir header [output_file_name] [reference]";
            print_header(argv[2],argv[3]);
        }
        else if (mode == "assemble") {
            if (argc != 9) throw "Usage:10 parameters needed\tpamir assemble [partition-file] [reference] [range] [output-file-vcf] [max-len] [read-length] dir_prefix";
            assemble(argv[2], argv[3], argv[4], argv[5], atoi(argv[6]), atoi(argv[7]), argv[8]);
        }

        else if (mode == "index_log") {
            if (argc != 3) throw "Usage:\tpamir index_log [log-file]\n\tIndexing the log file of assemble subcommand for fast random access. Required for running subcommand output_log.";
            log_idx( argv[2]);
        }
        else if (mode == "output_log") {
            if (argc != 4) throw "Usage:\tpamir log_idx [log-file] [range]\n\tRandom access to the log file of assemble subcommand. The index file generated by index_log subcommand is required, and the range argument x-y corresponds to [x,y). For example, range 10-20 corresponds to clusters from 10 to 19.";
            output_log( argv[2], argv[3]);
        }else if(mode == "process_range"){
            return range_process::main(argc-1,argv+1);
        }
        else if(mode == "process_orphan"){
            return orphan_process::main(argc-1,argv+1);
        }
        else if(mode == "recalibrate"){
            return recalibrate::main(argc-1,argv+1);
        }
        else if(mode == "smoother"){
            return smoother::main(argc-1,argv+1);
        }
        else if (mode == "sketch") {
            if (argc != 6) throw "Usage:4 parameters needed\tpamir sketch [long-read-file] [dat-path] [kmer-size] [window-size]";
            sketch(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]));
        }
        else if (mode == "find") {
            Logger::instance().debug.set_file("debug.find.log");
            if (argc != 7) throw "Usage:5 parameters needed\tpamir find [partition-file] [dat-files-path] [range] [reference] [short-read-len]";
            find_reads(argv[2], argv[3], argv[4], argv[5], stoi(argv[6]));
        }
        else if (mode == "extract") {
            Logger::instance().debug.set_file("debug.extract.log");
            if (argc != 8) throw "Usage:6 parameters needed\tpamir extract [partition-file] [long-read-file] [range] [p3-name] [dat_path] [ranges_file]";
            extract_reads(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
        }
        else if (mode == "consensus") {
            Logger::instance().debug.set_file("debug.consensus.log");
            if (argc != 10) throw "Usage:8 parameters needed\tpamir consensus [partition-file] [reference] [long-read-file] [dat-path] [range] [output-file-vcf] [max-len] dir_prefix";
            consensus(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], atoi(argv[8]), argv[9]);
        }
        else if (mode == "partition-hybrid") {
            if (argc != 4) throw "Usage:2 parameters needed\tpamir partition-hybrid [short-read-mapping] [short-read-len]";
            genome_partition_hybrid pt("partition", true);
            pt.cluster_reads(argv[2], stoi(argv[3]));
        }
        else {
            throw "Invalid mode selected";
        }
    }
    catch (const char *e) {
        Logger::instance().error("Error: %s\n", e);
        exit(1);
    }


    return 0;
}