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
/******************************************************************/
void append_vcf(const string &chrName, const string &reference, const vector< tuple< string, int, int, string, int, float > > &reports, const int &clusterId, string &vcf_str, string &vcf_str_del )
{
	for(int r=0;r<reports.size();r++)
	{
    	if(get<0>(reports[r])== "INS")
		{
			vcf_str += 	chrName;	vcf_str += 	"\t";
			vcf_str +=	std::to_string(get<1>(reports[r]));	vcf_str += "\t.\t";
			vcf_str +=  reference.at(r);
			vcf_str +=  "\t";
            vcf_str += reference.at(r);
            vcf_str += get<3>(reports[r]);
            vcf_str +=  "\t";
			vcf_str +=  std::to_string( get<5>(reports[r]));
			vcf_str +=  "\tPASS\t";
			vcf_str += "Cluster=";	vcf_str +=	std::to_string( clusterId ) ;
			vcf_str += ";Support=";	vcf_str	+=	std::to_string( get<4>(reports[r])) ;
			vcf_str += "\n";
		}
	}
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
void find_reads (const string &partition_file, const string &dat_path, const string &range)
{
    Sketch lr_sketch = Sketch(dat_path);

    //TODO FIX NAME - Get ouput as an argument
    string p2 = "partition-p2-" + range;
    genome_partition pt(partition_file, range);
    p2_partition pt_2(p2, true);

    while (1)
    {
        auto p 			= pt.read_partition();
		
        // end of the partition file
        if ( !p.size() )
            break;

        vector<pair<string, pair<pair<int, int>, pair<int, int> > > > cut_candidates;

        string chrName  = pt.get_reference();
        int pt_start    = pt.get_start();
        int pt_end      = pt.get_end();

        // cluster has too many or too few reads
        if ( p.size() > 7000 || p.size() <= 2 ) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID      : %d\n", pt.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.size());
            Logger::instance().info("INFO: Skipped Processing - Too few or Too many short reads\n");

            pt_2.add_cuts(p, cut_candidates, pt_start, pt_end, chrName);
            continue;
        }

        vector<string> reads;
        for (int i =0;i<p.size();i++)
        {
            reads.push_back(p[i].first.second);
        }

        lr_sketch.sketch_query(reads);

        cut_candidates = lr_sketch.find_cuts(true);

        if (cut_candidates.size() == 0) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID      : %d\n", pt.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.size());
            Logger::instance().info("INFO: Skipped Processing - No Long reads found\n");

            pt_2.add_cuts(p, cut_candidates, pt_start, pt_end, chrName);
            continue;
        }

        pt_2.add_cuts(p, cut_candidates, pt_start, pt_end, chrName);
    }
}
/*********************************************************************************************/
void extract_reads(const string &partition_file, const string &longread, const string &range, const string &p3_name)
{
	p2_partition pt(partition_file, range);
	cut_ranges ranges = cut_ranges(longread);
	
	while (1) {
		auto p 			= pt.read_partition();
	
		if ( !p.first.size() )
			break;

        // cluster has too many or too few reads
        if ( p.second.size() == 0 ) {
            continue;
        }

		for (int i = 0; i < p.second.size(); i++) {
            ranges.add_range(p.second[i].first, p.second[i].second.first.first, p.second[i].second.first.second);
		}
	}

    ranges.extract();
	
	//TODO FIX NAME + get output name as argument
	string p3 = "partition-p3-" + p3_name;
	p2_partition pt_2(partition_file, range);
	p3_partition pt_3(p3, true);
	
	while (1) {

		auto p 			= pt_2.read_partition();
		if ( !p.first.size() )
			break;

		vector<pair<pair<string, string>, pair<pair<int, int>, pair<int, int> > > > cuts;

        string chrName  = pt_2.get_reference();
        int pt_start    = pt_2.get_start();
        int pt_end      = pt_2.get_end();

        // cluster has too many or too few reads
        if ( p.first.size() > 7000 || p.first.size() <= 2 ) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID      : %d\n", pt_2.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.first.size());
            Logger::instance().info("INFO: Skipped Processing - Too few or Too many short reads\n");

            pt_3.add_reads(p.first, cuts, pt_start, pt_end, chrName);
            continue;
        }
        // cluster has too many or too few reads
        else if ( p.second.size() == 0 ) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID      : %d\n", pt_2.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.first.size());
            Logger::instance().info("INFO: Skipped Processing - No long reads found\n");

            pt_3.add_reads(p.first, cuts, pt_start, pt_end, chrName);
            continue;
        }

		for (int i = 0; i < p.second.size(); i++) {
            string tmp = ranges.get_cut(p.second[i].first, p.second[i].second.first.first,
                                        p.second[i].second.first.second);
            string cut;
            if (p.second[i].second.second.second == REV)
                cut = reverse_complement(tmp);
            else
                cut = tmp;
            cuts.push_back({{p.second[i].first, cut}, {{p.second[i].second.first.first, p.second[i].second.first.second},
                                                {p.second[i].second.second.first, p.second[i].second.second.second}}});
        }

		fix_reverse(cuts);
        pt_3.add_reads(p.first, cuts, pt_start, pt_end, chrName);
	}
}
/*********************************************************************************************/
void consensus (const string &partition_file, const string &reference, const string lr_path, const string dat_path, const string &range, const string &name, int max_len, const string &prefix)
{
	const double MAX_AT_GC 		= 0.7;
	const int MAX_REF_LEN		= 300000000;
	int LENFLAG					= 1000;//500;//1000;
	char *line 					= new char[MAX_CHAR];
    string out_vcf = prefix + "/" + name + ".vcf";
	FILE *fo_vcf 				= fopen(out_vcf.c_str(), "w");
    string out_vcf_del = prefix + "/" + name + "_DELS.vcf";
    string out_vcf_lq = prefix  + "/" + name + "_LOW_QUAL.vcf";
    FILE *fo_vcf_del 			= fopen(out_vcf_del.c_str(), "w");
	FILE *fo_vcf_lq 			= fopen(out_vcf_lq.c_str(), "w");
	
	genome ref(reference.c_str());
	map<string,string> chroms;
	p3_partition pt(partition_file, range);
	aligner al(max_len + 2010);
	auto alignment_engine = spoa::AlignmentEngine::Create(spoa::AlignmentType::kSW, 2, -32, -64, -1);
	spoa::Graph graph{};
	InsertionAssembler ia = InsertionAssembler(dat_path, lr_path);
	
	string tmp_ref; tmp_ref.reserve(4);
	string tmp_ref_lq; tmp_ref_lq.reserve(4);
	string vcf_info=""; vcf_info.reserve(10000000);
	string vcf_info_lq=""; vcf_info_lq.reserve(10000000);
	string vcf_info_del=""; vcf_info.reserve(10000000);
	const int MAX_BUFFER = 500;
	int n_buffer         =   0;
	int n_buffer2         =   0;

	int cnt = 0;
	char comment[20];
	ProgressBar progress(80);

	while (1) 
	{
		auto p 			= pt.read_partition();

		// end of the partition file
		if ( !p.first.size() )
			break;
	
        cnt += 1;
       	sprintf(comment, "%10d / %-10d", cnt, pt.get_total());
       	progress.update(((float)cnt/(float)pt.get_total()) * 100, comment);

		// cluster has too many or too few reads
		if ( p.first.size() > 7000 || p.first.size() <= 2 ) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID            : %d\n", pt.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.first.size());
            Logger::instance().info("INFO: Skipped Processing - Too few or Too many short reads\n");
            continue;
        }
		else if (p.second.size == 0) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID            : %d\n", pt.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.first.size());
            Logger::instance().info("INFO: Skipped Processing - No long reads found\n");
            continue;
		}
		else if ( p.second.size > 150 || p.second.size <= 2 ) {
            Logger::instance().info("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
            Logger::instance().info(" + Cluster ID            : %d\n", pt.get_id());
            Logger::instance().info(" + Short Reads Count     : %lu\n", p.first.size());
            Logger::instance().info(" + Long Reads Count      : %lu\n", p.second.size);
            Logger::instance().info("INFO: Skipped Processing - Too few or Too many long reads\n");
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
        Logger::instance().info(" + Cluster ID            : %d\n", cluster_id);
        Logger::instance().info(" + Short Reads Count     : %lu\n", p.first.size());
        Logger::instance().info(" + Long Reads Count      : %lu\n", p.second.size);
        Logger::instance().info(" + Spanning Range        : %s:%d-%d\n", chrName.c_str(), pt_start, pt_end);
        Logger::instance().info(" + Discovery Range       : %s:%d-%d\n", chrName.c_str(), ref_start, ref_end);
        Logger::instance().info(" + Reference             : %s\n", ref_part.c_str());
		
		// if the genomic region is too big
		if (ref_end - ref_start > MAX_REF_LEN) 
			continue;
		
		// holding the calls info, can be used to detect the repeated calls, etc.
		vector< tuple< string, int, int, string, int, float > > reports;//reports.clear();
		vector< tuple< string, int, int, string, int, float > > reports_lq;//reports.clear();

        vector<pair<string, int> > consensus;

		//BIMODAL
        if (p.second.bimodal) {
            Logger::instance().info(" + Type                  : bimodal\tleft: %d right: %d bimodal: %d\n\n", p.second.left_cuts.size(),
                                    p.second.right_cuts.size(), p.second.bimodal_cuts.size());

            graph.Clear();
            for (int i = 0; i < p.second.bimodal_cuts.size(); i++) {
                auto alignment = alignment_engine->Align(p.second.bimodal_cuts[i].first.second, graph);
                graph.AddAlignment(alignment, p.second.bimodal_cuts[i].first.second);
            }
            for (int i = 0; i < p.second.left_cuts.size(); i++) {
                auto alignment = alignment_engine->Align(p.second.left_cuts[i].first.second, graph);
                graph.AddAlignment(alignment, p.second.left_cuts[i].first.second);
            }
            for (int i = 0; i < p.second.right_cuts.size(); i++) {
                auto alignment = alignment_engine->Align(p.second.right_cuts[i].first.second, graph);
                graph.AddAlignment(alignment, p.second.right_cuts[i].first.second);
            }

            auto msa = graph.GenerateMultipleSequenceAlignment(true);

            pair<string, pair<int, int>> ans = cut_consensus_bimodal(msa, p.second.left_cuts.size(), p.second.right_cuts.size(),
                                                                     p.second.bimodal_cuts.size());

            consensus.push_back({ans.first, p.second.left_cuts.size() + p.second.bimodal_cuts.size() + p.second.right_cuts.size()});
        }
		else {
            //SINGLE PEAK
            if (p.second.right_cuts.size() < 2) {
                Logger::instance().info(" + Type              : single peak\n\n");
                
                graph.Clear();
                for (int i = 0; i < p.second.left_cuts.size(); i++) {
                    auto alignment = alignment_engine->Align(p.second.left_cuts[i].first.second, graph);
                    graph.AddAlignment(alignment, p.second.left_cuts[i].first.second);
                }

                auto msa = graph.GenerateMultipleSequenceAlignment(true);

                pair<string, pair<int, int>> ans = cut_consensus_single(msa);

                consensus.push_back({ans.first, p.second.left_cuts.size()});
            }
            //LONG INSERTION
            else if (!p.second.left_cuts.empty() && !p.second.right_cuts.empty()) {
                Logger::instance().info(" + Type              : long insertion\n\n");

                vector<string> reads_1, reads_2;
                vector<string> left_reads, right_reads;

                graph.Clear();
                for (int i = 0; i < p.second.left_cuts.size(); i++) {
                    auto alignment = alignment_engine->Align(p.second.left_cuts[i].first.second, graph);
                    graph.AddAlignment(alignment, p.second.left_cuts[i].first.second);
                    reads_1.push_back(p.second.left_cuts[i].first.second);
                }
                auto msa_1 = graph.GenerateMultipleSequenceAlignment(true);
                string cons_1 = msa_1[msa_1.size() - 1];
                cons_1.erase(std::remove(cons_1.begin(), cons_1.end(), '-'), cons_1.end());

                graph.Clear();
                for (int i = 0; i < p.second.right_cuts.size(); i++) {
                    auto alignment = alignment_engine->Align(p.second.right_cuts[i].first.second, graph);
                    graph.AddAlignment(alignment, p.second.right_cuts[i].first.second);
                    reads_2.push_back(p.second.right_cuts[i].first.second);
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
    
                pair<string, int> ans = ia.assemble(left_reads, right_reads);

                consensus.push_back(ans);
            }
        }

        for (int i = 0; i < consensus.size(); i++) {
            int contig_support = consensus[i].second;

            Logger::instance().info("\n\n>>>>> Length: %d Support: %d Contig: %s\n", consensus[i].first.size(),
                                    contig_support, consensus[i].first.c_str());

            //Align consensus
            al.align(ref_part, consensus[i].first);
            if(al.extract_calls(cluster_id, reports_lq, reports, contig_support, ref_start,">>>")==0) {
                string rc_contig = reverse_complement(consensus[i].first);
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
            if (argc != 5) throw "Usage:3 parameters needed\tpamir find [partition-file] [dat-files-path] [range]";
            find_reads(argv[2], argv[3], argv[4]);
        }
		else if (mode == "extract") {
			if (argc != 6) throw "Usage:4 parameters needed\tpamir extract [partition-file] [long-read-file] [range] [p3-name]";
			extract_reads(argv[2], argv[3], argv[4], argv[5]);
		}
		else if (mode == "consensus") {
            if (argc != 10) throw "Usage:8 parameters needed\tpamir consensus [partition-file] [reference] [long-read-file] [dat-path] [range] [output-file-vcf] [max-len] dir_prefix";
            consensus(argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], atoi(argv[8]), argv[9]);
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