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
#include "assembler_ext.h"
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
void modifyOeaUnmapped (const string &path, const string &result)
{
	ifstream fin(path.c_str());
	FILE *fo = fopen(result.c_str(), "w");
	
	string l, a, b, aa, bb, prevaa; 
	while (getline(fin, a)) {
		getline(fin, b);
		getline(fin, l);
		getline(fin, l);

		prevaa=aa;
		int ll = a.size();
		if (ll > 1 && a[ll - 2] == '/')
			a = a.substr(0, ll - 2);
		if (a == aa)
			continue;
		else {
			if(aa!="")
				fprintf(fo, "%s %s\n", aa.c_str(), bb.c_str());
			aa = a, bb = b;
		}
	}
	if(prevaa != aa)
	{
		fprintf(fo, "%s %s\n", aa.c_str(), bb.c_str());
	}
	fin.close();
	fclose(fo);
}
/****************************************************************/
void partify (const string &read_file, const string &mate_file, const string &out, int threshold) 
{
	gzFile gzf = gzopen(mate_file.c_str(), "rb");
	unordered_map<string, string> mymap;
	const int MAXB = 259072;
	char buffer[MAXB];
	char name[MAXB], read[MAXB];
	while (gzgets(gzf, buffer, MAXB)) {
		sscanf(buffer, "%s %s", name, read);
		string read_name = read;
		if (read_name.size() > 2 && read_name[read_name.size() - 2] == '/')
			read_name = read_name.substr(0, read_name.size() - 2);
		mymap[string(name + 1)] = read_name;
	}
	auto g=mymap.begin();
	while(g!=mymap.end())
	{
		g++;
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
	fclose(fo);
	fclose(fidx);
}
/****************************************************************/
// For outputing specific log
void log_idx (const string &log_file ) 
{
	FILE *fin = fopen(log_file.c_str(), "rb");
	FILE *fidx = fopen((log_file + ".idx").c_str(), "wb");
	char *readline = (char*)malloc(50000);
	char *token = (char*)malloc(100);
	size_t idx_pos = ftell(fin);
	int l_id, offset;	
	int num_inserted = 0; // to resolve skipping partition issue

	fwrite( &idx_pos, 1, sizeof(size_t), fidx); // initialize an log for partition id ZERO
	while( NULL != fgets( readline, 50000, fin ) )
	{	
		if ( 0 == strncmp("PARTITION ", readline, 10) )
		{
			sscanf(readline, "%s %d %n", token, &l_id, &offset);
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
			end = tok ? atol(tok) : -1;
		}
		free(dup);
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
	char pref[MAXB];
	char name[MAXB], read[MAXB];
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
	if ( 0 != strncmp("PARTITION ", pref, 10) )
	{	exit(1); fprintf(stderr, "Incorrect Start at %s", pref);
	}
	fprintf( fo, "%s", pref);
	
	fgets(pref, MAXB, fi);
	while ( 0 != strncmp("PARTITION ", pref, 10) )
	{
		fprintf( fo, "%s", pref);
		fgets(pref, MAXB, fi);
	}

	num_read = 0;

	fclose(fi);
	fclose(fo);
	if ( num_read == 0)
		goto reset;
	return num_cluster;
}
/******************************************************************/
string assemble_with_sga (const string &input_fastq)
{
	char *sgapython = new char[100000];
	strcpy(sgapython,"/cs/compbio3/yenyil/Pinar/pinsertionForExperiments/sga.py");
	char *sgacmd = new char[100000];
	sprintf(sgacmd, "python %s %s",sgapython,input_fastq.c_str());
	system(sgacmd);
	string outofsga = input_fastq+string(".sgaout.fa");
	return outofsga;
}
/********************************************************************************/
string prepare_sga_input (const string &out_vcf, const vector<pair<pair<string, string>, pair<int, int>>> &p, const int &read_length)
{
	string qual(read_length,'I');
	string inputforsga = string(out_vcf + "_fastq.fq");
	FILE *fqforsga 	   = fopen(inputforsga.c_str(),"w");
	for(int i =0;i < p.size(); i++)
	{
		if(p[i].first.second.length()>read_length)
		{
			int j		= 1;
			int stpoint = 0;
			while( stpoint + read_length < p[i].first.second.length() )
			{
				string read = p[i].first.second.substr( stpoint, read_length );
				fprintf( fqforsga,"@%s\n%s\n+\n%s\n",( p[i].first.first + "_" + itoa(j) ).c_str(), read.c_str(), qual.c_str() );
				j++;
				stpoint += 10;
			}
			string read = p[i].first.second.substr( p[i].first.second.length() - read_length, read_length );
			fprintf( fqforsga,"@%s\n%s\n+\n%s\n",( p[i].first.first + "_" + itoa(j) ).c_str(), read.c_str(), qual.c_str() );
		}
		else
			fprintf( fqforsga, "@%s\n%s\n+\n%s\n", p[i].first.first.c_str(), p[i].first.second.c_str(), qual.c_str() );
	}
	fclose(fqforsga);
	return inputforsga;
}
/*****************************************************************/
void print_calls(string chrName, vector< tuple< string, int, int, string, int, float > > &reports, FILE *fo_vcf, const int &clusterId)
{
	for(int r=0;r<reports.size();r++)
	{
		if(get<0>(reports[r])== "INS")
			fprintf(fo_vcf, "%s\t%d\t%d\t.\t%s\t%f\tPASS\tSVTYPE=INS;LEN=%d;END=%d;Cluster=%d;Support=%d;Identity=%f\n", chrName.c_str(), get<1>(reports[r]), get<2>(reports[r]), get<3>(reports[r]).c_str(), -10*log(1-get<5>(reports[r])), get<2>(reports[r]), get<1>(reports[r]), clusterId, get<4>(reports[r]), get<5>(reports[r]) ); 
	}
}
/****************************************************************/
void assemble (const string &partition_file, const string &reference, const string &range, const string &out_vcf, int max_len, int read_length, const int &hybrid)
{
	const double MAX_AT_GC 		= 0.7;
	const int ANCHOR_SIZE 		= 16;
	const int MAX_REF_LEN		= 300000000;
	int LENFLAG					= 1000;
	char *line 			= new char[1000000];
	FILE *fo_vcf 				= fopen(out_vcf.c_str(), "w");
	
	assembler as(max_len, 15);
	genome ref(reference.c_str());
	genome_partition pt;
	
	aligner al( ANCHOR_SIZE, max_len + 2010 );
	while (1) 
	{
		auto p 			= pt.read_partition(partition_file, range);
		// end of the partition file
		if ( !p.size() ) 
			break;
		
		// cluster has too many or too few reads
		if ( p.size() > 100000 || p.size() <= 2 ) 	
			continue;
		
		string chrName  = pt.get_reference();
		int cluster_id  = pt.get_cluster_id();
		int pt_start    = pt.get_start();
		int pt_end      = pt.get_end();
		int ref_start   = pt_start - LENFLAG;
		int ref_end     = pt_end   + LENFLAG;
		string ref_part = ref.extract(chrName, ref_start, ref_end);
		fprintf(stdout,"-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-*-<=*=>-\n");
		fprintf(stdout," + Cluster ID      : %d\n", cluster_id);
		fprintf(stdout," + Reads Count     : %lu\n", p.size());
		fprintf(stdout," + Spanning Range  : %s:%d-%d\n", chrName.c_str(), pt_start, pt_end);
		fprintf(stdout," + Discovery Range : %s:%d-%d\n\n", chrName.c_str(), ref_start, ref_end);

		// if the genomic region is too big
		if (ref_end - ref_start > MAX_REF_LEN) 
			continue;

		// holding the calls info, can be used to detect the repeated calls, etc.
		vector< tuple< string, int, int, string, int, float > > reports;


		auto contigs    = as.assemble(p);
		for ( auto &contig: contigs )
		{
			int contig_support		= contig.support();
			int con_len 			= contig.data.length();
			if( check_AT_GC(contig.data, MAX_AT_GC) == 0 || contig_support <= 1 || con_len > max_len + 400 ) continue;
		
			fprintf(stdout, "\n>>>>> Length: %d Support: %d\n", con_len, contig_support);
			for(int z=0;z<contig.read_information.size();z++)
				fprintf(stdout,"%s %d %d %s\n", 
					contig.read_information[z].name.c_str(),
					contig.read_information[z].location,
					contig.read_information[z].in_genome_location, 
					contig.read_information[z].data.c_str());
			
			fprintf (stdout, "\n--->\n");
			al.align(ref_part, contig.data);
			if(al.extract_calls(cluster_id, reports, contig_support, ref_start)==0)
			{
				string rc_contig = reverse_complement(contig.data);	
				fprintf (stdout, "\n<---\n");
				al.align(ref_part, rc_contig);
				al.extract_calls(cluster_id, reports, contig_support, ref_start);
			}
		}
		print_calls(chrName, reports, fo_vcf, pt.get_cluster_id());
		if( ( reports.size() == 0 || reports.size() > 1 ) && hybrid == 1)
		{
			reports.clear();
			string outofsga 	= assemble_with_sga( prepare_sga_input( out_vcf, p, read_length ) );
			FILE *fcontig 		= fopen(outofsga.c_str(),"r");
			int contig_support 	= p.size();
			while( fgets( line, 1000000, fcontig ) != NULL )
			{
				fgets( line, 1000000, fcontig );
				line[ strlen(line)-1 ]		='\0';
				string contig 				= string(line);
				int con_len 				= contig.length();
				if( check_AT_GC( contig, MAX_AT_GC ) == 0 || contig_support <=1 || con_len > max_len + 400 ) continue;
				al.align(ref_part, contig);
				if(al.extract_calls(cluster_id, reports, contig_support, ref_start)==0)
				{
					string rc_contig = reverse_complement(contig);	
					al.align(ref_part, rc_contig);
					al.extract_calls(cluster_id, reports, contig_support, ref_start);
				}
			}
			print_calls( chrName, reports, fo_vcf, pt.get_cluster_id());
			reports.clear();
		}

	}
	fclose(fo_vcf);
}
/*********************************************************************************************/
int main(int argc, char **argv)
{
	try {
		if (argc < 2) throw "Usage:\tsniper [mode=(?)]";

		string mode = argv[1];
		if (mode == "verify_sam") {
			if (argc < 7) throw "Usage:\tsniper verify_sam [sam-file] [output] outputtype oea? orphan?";
			extractor ext(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
		}
		else if (mode == "remove_concordant") {
			if (argc != 7) throw "Usage:\tsniper remove_concordant [sam-file] [output] outputtype oea? orphan?";
			extractor ext(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
		}
		else if (mode == "mask" || mode == "maski") {
			if (argc != 6) throw "Usage:\tsniper mask/maski [repeat-file] [reference] [output] [padding]";
			mask(argv[2], argv[3], argv[4], atoi(argv[5]), mode == "maski");
		}
		else if (mode == "sort") {
			if (argc != 4) throw "Usage:\tsniper sort [sam-file] [output]";
			sortFile(argv[2], argv[3], 2 * GB);
		}
		else if (mode == "modify_oea_unmap") {
			if (argc != 4) throw "Usage:\tsniper modi_oea_unmap [fq-file] [output]";
			modifyOeaUnmapped(argv[2], argv[3]);
		}
		else if (mode == "partition") {
			if (argc != 6) throw "Usage:\tsniper partition [read-file] [mate-file] [output-file] [threshold]";
			partify(argv[2], argv[3], argv[4], atoi(argv[5]));
		}
		else if (mode == "assemble") {
			if (argc != 9) throw "Usage:10 parameters needed\tsniper assemble [partition-file] [reference] [range] [output-file-vcf] [max-len] [read-length] [hybrid]"; 
			assemble(argv[2], argv[3], argv[4], argv[5], atoi(argv[6]), atoi(argv[7]), atoi(argv[8]));
		}
		else if (mode == "get_cluster") {
			if (argc != 4) throw "Usage:\tsniper get_cluster [partition-file] [range]";
			genome_partition pt;
			pt.output_partition( argv[2], argv[3]);
		}
		else if (mode == "log_idx") {
			if (argc != 3) throw "Usage:\tsniper log_idx [log-file]";
			log_idx( argv[2]);
		}
		else if (mode == "output_log") {
			if (argc != 4) throw "Usage:\tsniper log_idx [log-file] [range]";
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
