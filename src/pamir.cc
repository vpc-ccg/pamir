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
#include "bam_parser.h"
#include "sort.h"


using namespace std;

inline string space (int i) 
{
	return string(i, ' ');
}
/********************************************************************/
inline string itoa (int i)
{
	char c[50];
	sprintf(c, "%d", i);
	return string(c);
}
/*******************************************************************/
void mask (const string &repeats, const string &path, const string &result, int pad = 0, bool invert = false)
{
	const int MAX_BUF_SIZE = 2000;
	size_t b, e, m = 0;
	char ref[MAX_BUF_SIZE], name[MAX_BUF_SIZE], l[MAX_BUF_SIZE];
	
	map<string, vector<pair<size_t, size_t>>> masks;
	FILE *fi = fopen(repeats.c_str(), "r");
	int prev = 0; string prev_ref = "";
	while (fscanf(fi, "%s %lu %lu %s", ref, &b, &e, name) != EOF) {
		if (e - b < 2 * pad) continue;
		masks[ref].push_back({b + pad, e - pad});
	}
	for (auto &r: masks) 
		sort(r.second.begin(), r.second.end());
	if (invert) for (auto &r: masks) 
	{
		vector<pair<size_t, size_t>> v;
		size_t prev = 0;
		for (auto &p: r.second) {
			if (prev < p.first)
				v.push_back({prev, p.first});
			prev = p.second;
		}
		v.push_back({prev, 500 * 1024 * 1024});
		r.second = v;
	}
	fclose(fi);
	
	const int MAX_CHR_SIZE = 300000000;
	char *x = new char[MAX_CHR_SIZE];

	fi = fopen(path.c_str(), "r");
	FILE *fo = fopen(result.c_str(), "w");
	fgets(l, MAX_BUF_SIZE, fi);
	
	while (true) {
		if (feof(fi))
			break;
		
		fprintf(fo, "%s", l);
		char *p = l; while (*p && !isspace(*p)) p++; *p = 0;
		string n = string(l + 1);
		printf("Parsed %s\n", n.c_str());
		
		size_t xlen = 0;
		while (fgets(l, MAX_BUF_SIZE, fi)) {
			if (l[0] == '>') 
				break;
			if(l[strlen(l)-1]=='\n'){
				memcpy(x + xlen, l, strlen(l)- 1);
				xlen += strlen(l)-1;
			}
			else{
				memcpy(x + xlen, l, strlen(l));
				xlen += strlen(l);
			}
		}
		x[xlen] = 0;

		for (auto &r: masks[n]) {
			if (r.first >= xlen) continue;
			if (r.second >= xlen) r.second = xlen;
			memset(x + r.first, 'N', r.second - r.first);
		}

		for (size_t i = 0; i < xlen; i += 120) 
			fprintf(fo, "%.120s\n", x + i);
	}

	fclose(fi);
	fclose(fo);
	delete[] x;
}
/**************************************************************/
void partify (const string &read_file, const string &out, int threshold, const string &mate_file) 
{
	FILE *fin = fopen(mate_file.c_str(), "r");
	unordered_map<string, string> mymap;
	const int MAXB = 259072;
	char name[MAXB], read[MAXB], tmp[MAXB];
	while (fgets(name, MAXB, fin)) {
		fgets(read, MAXB, fin);
		fgets(tmp, MAXB, fin);
		fgets(tmp, MAXB, fin);
		if (strlen(name)>2 && name[strlen(name) - 3] == '/')
			name[strlen(name)-3]='\0';
		read[strlen(read)-1]='\0';
		mymap[string(name+1)] = read;
	}
	genome_partition pt(read_file, threshold, mymap);
	int fc = 1;
	FILE *fo = fopen(out.c_str(), "wb");
	FILE *fidx = fopen((out + ".idx").c_str(), "wb");
	while (pt.has_next()) {
		auto p = pt.get_next();
		size_t i = pt.dump(p, fo, fc);
		fwrite(&i, 1, sizeof(size_t), fidx);
		fc++;
	}
	fclose(fin);
	fclose(fo);
	fclose(fidx);
}
void partify_orphan (const string &read_file, const string &out, int threshold, 
	const string &contig_file, const string &oea2orphan,const string &mate_file) 
{
	FILE *fin = fopen(mate_file.c_str(), "r");
	unordered_map<string, string> mymap;
	const int MAXB = 259072;
	char name[MAXB], read[MAXB], tmp[MAXB];
	while (fgets(name, MAXB, fin)) {
		fgets(read, MAXB, fin);
		fgets(tmp, MAXB, fin);
		fgets(tmp, MAXB, fin);
		if (strlen(name)>2 && name[strlen(name) - 3] == '/')
			name[strlen(name)-3]='\0';
		read[strlen(read)-1]='\0';
		mymap[string(name+1)] = read;
	}
	genome_partition pt(read_file, threshold, mymap); 
	// loading orphan association information
	pt.load_orphan( contig_file, oea2orphan);
	
	int fc = 1;
	FILE *fo = fopen(out.c_str(), "wb");
	FILE *fidx = fopen((out + ".idx").c_str(), "wb");
	while (pt.has_next()) {
		auto p = pt.get_next();
		size_t i = pt.dump(p, fo, fc);
		fwrite(&i, 1, sizeof(size_t), fidx);
		fc++;
	}
	fprintf(stdout, "%d\n",fc-1);
	fclose(fin);
	fclose(fo);
	fclose(fidx);
}

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
void print_calls(string chrName, const string &reference, vector< tuple< string, int, int, string, int, float > > &reports, FILE *fo_vcf, const int &clusterId)
{
	for(int r=0;r<reports.size();r++)
	{
		if(get<0>(reports[r])== "INS"){
			fprintf(fo_vcf, "%s\t",	 			chrName.c_str());
			fprintf(fo_vcf, "%d\t", 			get<1>(reports[r]));
		//	fprintf(fo_vcf, "%d\t", 			get<2>(reports[r]));
			fprintf(fo_vcf, ".\t");
			fprintf(fo_vcf, "%c\t",reference.at(r));
			fprintf(fo_vcf, "<INS>\t");
		//	fprintf(fo_vcf, "%s\t",				get<3>(reports[r]).c_str());
		//	fprintf(fo_vcf, "%f\t", 			-10*log(1-get<5>(reports[r])));
			fprintf(fo_vcf, "%f\t", 			get<5>(reports[r]));
			fprintf(fo_vcf, "PASS\tSVTYPE=INS;");
			fprintf(fo_vcf,	"SVLEN=%d;",  		get<2>(reports[r]));
			fprintf(fo_vcf, "END=%d;",  		get<1>(reports[r]) + get<2>(reports[r])-1);
			fprintf(fo_vcf, "Cluster=%d;", 		clusterId);
			fprintf(fo_vcf, "Support=%d;", 		get<4>(reports[r]));
		//	fprintf(fo_vcf, "Identity=%f\t", 	get<5>(reports[r])); 
			fprintf(fo_vcf, "SEQ=%s\n", 		get<3>(reports[r]).c_str()); 
		}
	}
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
//			vcf_str +=  std::to_string( get<2>(reports[r])) ;
//			vcf_str += 	";END=";	vcf_str +=	std::to_string( get<1>(reports[r]) + get<2>(reports[r])-1 );
			vcf_str += "Cluster=";	vcf_str +=	std::to_string( clusterId ) ;
			vcf_str += ";Support=";	vcf_str	+=	std::to_string( get<4>(reports[r])) ;
//			vcf_str += ";SEQ=";		vcf_str	+=	get<3>(reports[r]);
			vcf_str += "\n";
		}
/*
		if(get<0>(reports[r])== "INS")
		{
			vcf_str += 	chrName;	vcf_str += 	"\t";
			vcf_str +=	std::to_string(get<1>(reports[r]));	vcf_str += "\t.\t";
			vcf_str +=  reference.at(r);
			vcf_str +=  "\t<INS>\t";
			vcf_str +=  std::to_string( get<5>(reports[r]));
			vcf_str +=  "\tPASS\tSVTYPE=INS;SVLEN=";
			vcf_str +=  std::to_string( get<2>(reports[r])) ;
			vcf_str += 	";END=";	vcf_str +=	std::to_string( get<1>(reports[r]) + get<2>(reports[r])-1 );
			vcf_str += ";Cluster=";	vcf_str +=	std::to_string( clusterId ) ;
			vcf_str += ";Support=";	vcf_str	+=	std::to_string( get<4>(reports[r])) ;
			vcf_str += ";SEQ=";		vcf_str	+=	get<3>(reports[r]);
			vcf_str += "\n";
		}
    	if(get<0>(reports[r])== "DEL")
		{
			vcf_str_del += 	chrName;	vcf_str_del += 	"\t";
			vcf_str_del +=	std::to_string(get<1>(reports[r]));	vcf_str_del += "\t.\t";
			vcf_str_del +=  reference.at(r);
			vcf_str_del +=  "\t<DEL>\t";
			vcf_str_del +=  std::to_string( get<5>(reports[r]));
			vcf_str_del +=  "\tPASS\tSVTYPE=DEL;SVLEN=";
			vcf_str_del +=  std::to_string( get<2>(reports[r])) ;
			vcf_str_del += 	";END=";	vcf_str_del +=	std::to_string( get<1>(reports[r]) + get<2>(reports[r])-1 );
			vcf_str_del += ";Cluster=";	vcf_str_del +=	std::to_string( clusterId ) ;
			vcf_str_del += ";Support=";	vcf_str_del	+=	std::to_string( get<4>(reports[r])) ;
			vcf_str_del += ";SEQ=";		vcf_str_del	+=	get<3>(reports[r]);
			vcf_str_del += "\n";
		}*/
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
	genome_partition pt;
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
		auto p 			= pt.read_partition(partition_file, range); 
		// end of the partition file
		if ( !p.size() ) 
			break;
		
		// cluster has too many or too few reads
		if ( p.size() > 7000 || p.size() <= 2 ) 
			continue;
		string chrName  = pt.get_reference();
		int cluster_id  = pt.get_cluster_id();
		int pt_start    = pt.get_start();
		int pt_end      = pt.get_end();
		int ref_start   = pt_start - LENFLAG;
		int ref_end     = pt_end   + LENFLAG;
		string ref_part = ref.extract(chrName, ref_start, ref_end);
		log("-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
		log(" + Cluster ID      : %d\n", cluster_id);
		log(" + Reads Count     : %lu\n", p.size());
		log(" + Spanning Range  : %s:%d-%d\n", chrName.c_str(), pt_start, pt_end);
		log(" + Discovery Range : %s:%d-%d\n", chrName.c_str(), ref_start, ref_end);
		log(" + Reference       : %s\n\n", ref_part.c_str());
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
		
			log("\n\n>>>>> Length: %d Support: %d Contig: %s\n", con_len, contig_support, contig.data.c_str());
			for(int z=0;z<contig.read_information.size();z++)
				log("%s %s %d %d\n", 
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
			tmp_ref += ref.getchar(chrName, tmp_end);
			//tmp_ref += ref.extract(chrName, tmp_end, tmp_end);
		}
		
		tmp_ref_lq.clear();//string tmp_ref = ""; 
		for (int j =0; j <reports_lq.size();j++)
		{
			int tmp_end = get<1>(reports_lq[j]);
			tmp_ref_lq += ref.getchar(chrName, tmp_end);
		}
		append_vcf( chrName, tmp_ref, reports, pt.get_cluster_id(), vcf_info, vcf_info_del);
		n_buffer++;
		if ( 0 == n_buffer%MAX_BUFFER )
		{
			fprintf( fo_vcf, "%s", vcf_info.c_str());
			n_buffer = 0;
			vcf_info.clear();
		}
		append_vcf( chrName, tmp_ref_lq, reports_lq, pt.get_cluster_id(), vcf_info_lq, vcf_info_del);
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
	string log_path = "";


	try {

		if (argc < 2) throw "Usage:\tpamir [mode=(?)]";

		string mode = argv[1];
		if (mode == "verify_sam") {
			if ( 4 == argc )
			{
				extractor ext(argv[2], argv[3], 0.95);
			}
			else if ( 5 == argc)
			{
				extractor ext(argv[2], argv[3], stod(argv[4]) );
			}
			else
			{	
				throw "Usage:\tpamir verify_sam [sam-file] [output_prefix] matched_ratio";
			}
		}
		else if (mode == "remove_concordant") {
			if (argc != 8) throw "Usage:\tpamir remove_concordant [sam-file] [output] outputtype oea? orphan? matched_ratio";
			extractor ext(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), stod(argv[7]));
		}
		else if (mode == "getfastq") {
			if(argc != 4) throw "Usage:\tpamir getfastq  [sam-file/bam-file] [output]";
			extractor ext(argv[2], argv[3] );//, atoi(argv[4]), atoi(argv[5]), atoi(argv[6]), stod(argv[7]));
		}
		else if (mode == "mask" || mode == "maski") {
			if (argc != 6) throw "Usage:\tpamir mask/maski [repeat-file] [reference] [output] [padding]";
			mask(argv[2], argv[3], argv[4], atoi(argv[5]), mode == "maski");
		}
		else if (mode == "sort") {
			if (argc != 4) throw "Usage:\tpamir sort [sam-file] [output]";
			sortFile(argv[2], argv[3], 2 * GB);
		}
		else if (mode == "partition") {
			if ( 6 == argc )
            {

				log_path = argv[5]; log_path += "_partitionprocessor.log";

				log_init( log_path );
                partify(argv[2], argv[3], atoi(argv[4]), argv[5]);


			}
			else if ( 8 == argc)
			{
				log_path = argv[5]; log_path += "_partitionprocessor.log";
				log_init( log_path );
				partify_orphan(argv[2], argv[3], atoi(argv[4]), argv[5], argv[6], argv[7]);
				log_close();
			}
			else{ throw "Usage:\tpamir partition [read-file] [output-file] [threshold] [ [orphan-contig] [oea2orphan] ] [mate_file]"; }
		}
		else if (mode=="header"){
			if (argc !=4) throw "Usage:3 parameters needed\tpamir header [output_file_name] [reference]";
			print_header(argv[2],argv[3]);
		}
		else if (mode == "assemble") {
			if (argc != 9) throw "Usage:10 parameters needed\tpamir assemble [partition-file] [reference] [range] [output-file-vcf] [max-len] [read-length] dir_prefix";
			log_path = argv[5]; log_path += ".log";
			log_init( "" );	//log_init( log_path );
			assemble(argv[2], argv[3], argv[4], argv[5], atoi(argv[6]), atoi(argv[7]), argv[8]);
			log_close();
		}
		else if (mode == "get_cluster") {
			if (argc != 4) throw "Usage:\tpamir get_cluster [partition-file] [range]";
			genome_partition pt;
			pt.output_partition( argv[2], argv[3]);
		}
		else if (mode == "index_log") {
			if (argc != 3) throw "Usage:\tpamir index_log [log-file]";
			log_idx( argv[2]);
		}
		else if (mode == "output_log") {
			if (argc != 4) throw "Usage:\tpamir log_idx [log-file] [range]";
			output_log( argv[2], argv[3]);
		}
		else {
			throw "Invalid mode selected";
		}
	}
	catch (const char *e) {
		ERROR("Error: %s\n", e);
		exit(1);
	}

		
	return 0;
}
