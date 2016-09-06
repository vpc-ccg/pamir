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
/*****************************************************************/
inline string itoa (int i)
{
	char c[50];
	sprintf(c, "%d", i);
	return string(c);
}
/*******************************************************************/
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
/********************************************************************************/
string assemble_with_sga (const string &out_vcf, const vector<pair<pair<string, string>, pair<int, int>>> &p, const int &read_length)
{
	string qual(read_length,'I');
	char *inputforsga = new char[10000];
	strcpy(inputforsga,(out_vcf+"_fastq.fq").c_str());
	FILE *fqforsga = fopen(inputforsga,"w");
	for(int i =0;i<p.size();i++)
	{
		if(p[i].first.second.length()>read_length)
		{
			int j=1;
			int stpoint=0;
			while(stpoint+read_length<p[i].first.second.length())
			{
				string read = p[i].first.second.substr(stpoint,read_length);
				fprintf(fqforsga,"@%s\n%s\n+\n%s\n",(p[i].first.first + "_" + itoa(j)).c_str(), read.c_str(),qual.c_str());
				j++;
				stpoint+=10;
			}
			string read = p[i].first.second.substr(p[i].first.second.length()-read_length,read_length);
			fprintf(fqforsga,"@%s\n%s\n+\n%s\n",(p[i].first.first + "_" + itoa(j)).c_str(), read.c_str(),qual.c_str());
		}
		else
			fprintf(fqforsga,"@%s\n%s\n+\n%s\n",p[i].first.first.c_str(),p[i].first.second.c_str(),qual.c_str());
	}
	fclose(fqforsga);
	char *sgapython = new char[100000];
	strcpy(sgapython,"/cs/compbio3/yenyil/Pinar/pinsertionForExperiments/sga.py");
	char *sgacmd = new char[100000];
	sprintf(sgacmd, "python %s %s",sgapython,inputforsga);
	system(sgacmd);
	string outofsga = string(inputforsga)+string(".sgaout.fa");
	return outofsga;
}
/******************************************************************/
string assemble_orphan_with_sga (const string &orphan_fastq)
{
	char *sgapython = new char[100000];
	strcpy(sgapython,"/cs/compbio3/yenyil/Pinar/pinsertionForExperiments/sga.py");
	char *sgacmd = new char[100000];
	sprintf(sgacmd, "python %s %s",sgapython,orphan_fastq.c_str());
	system(sgacmd);
	string outofsga = orphan_fastq+string(".sgaout.fa");
	return outofsga;
}
/*****************************************************************/
void print_calls(vector< tuple< string, int, int, string, int, float > > &reports, FILE *fo_vcf, FILE *fo_full, const int &clusterId)
{
	for(int r=0;r<reports.size();r++)
	{
		fprintf(fo_full,"%s\t%d\t%d\t%s\t%d\n------------------------------------------------------\n",get<0>(reports[r]).c_str(), get<1>(reports[r]), get<2>(reports[r]), get<3>(reports[r]).c_str(),get<4>(reports[r]));
//		fprintf(fo_vcf,"%s\t%d\t%d\t%s\t%d\n", get<0>(reports[r]).c_str(), get<1>(reports[r]), get<2>(reports[r]), get<3>(reports[r]).c_str(), get<4>(reports[r]));					
		fprintf(fo_vcf, "%s\t%d\t.\t%s\t%f\tPASS\tSVTYPE=INS;LEN=%d;END=%d;Cluster=%d;Support=%d;Identity=%f\n", get<0>(reports[r]).c_str(), get<1>(reports[r]), get<3>(reports[r]).c_str(), -10*log(1-get<5>(reports[r])), get<2>(reports[r]), get<1>(reports[r]), clusterId, get<4>(reports[r]), get<5>(reports[r]) ); 
	}
}
/****************************************************************/
void assemble (const string &partition_file, const string &reference, const string &range, const string &out_vcf, const string &out_full, int max_len, int read_length, const int &hybrid)
{
	const double MAX_AT_GC = 0.7;
	const int ANCHOR_SIZE = 16;
	const int MAX_REF_LEN=300000000;
	int LENFLAG=1000;
	int pt_start=-1;
	int pt_end=-1;
	FILE *fo_vcf = fopen(out_vcf.c_str(), "w");
	FILE *fo_vcf_del = fopen((out_vcf+"_del").c_str(), "w");
	FILE *fo_full = fopen(out_full.c_str(), "w");	
	FILE *fo_full_del = fopen((out_full+".del").c_str(), "w");	
	
	assembler as(max_len, 15);
	genome ref(reference.c_str());
	genome_partition pt;
	while (1) 
	{
		auto p = pt.read_partition(partition_file, range);
		if (!p.size()) 
			break;
		int cluster_id = pt.get_cluster_id();
		if (p.size() > 100000 || p.size() <= 2) 	continue;
		pt_start=pt.get_start();
		pt_end=pt.get_end();
		auto contigs = as.assemble(p);
		fprintf(fo_full,"PARTITION %d | read number in partition= %lu; contig number in partition= %lu; partition range=%d\t%d\n*******************************************************************************************\n",cluster_id, p.size(), contigs.size(), pt_start, pt_end);
		int contigNum=1;
		vector< tuple< string, int, int, string, int, float > > reports;
		int contigSupport;
		for (auto &contig: contigs)
		{
			contigSupport=contig.support();
			int con_len = contig.data.length();
			if(check_AT_GC(contig.data, MAX_AT_GC)==0 || contigSupport <=1 || con_len > max_len + 400 || (pt_end +1000 - (pt_start - 1000)) > MAX_REF_LEN ) continue;
			
			fprintf(fo_full,"#################################################################################################\nNEWCONTIG %d IN PARTITION %d length %lu; readNum %d; contig range %d\t%d\n%s\nSUPPORTIVE READS ARE:\n", contigNum, cluster_id, contig.data.length(), contigSupport, pt_start, pt_end, contig.data.c_str());
			for(int z=0;z<contig.read_information.size();z++)
				fprintf(fo_full,"%s %d %d %s\n",contig.read_information[z].name.c_str(),contig.read_information[z].location,contig.read_information[z].in_genome_location, contig.read_information[z].data.c_str());
			string ref_part = ref.extract(pt.get_reference(), pt_start - LENFLAG, pt_end + LENFLAG);
			int matrix_size = ref_part.size();
			if(con_len > matrix_size) matrix_size = con_len;
			aligner al(ANCHOR_SIZE, matrix_size);
			fprintf(fo_full, "REF region is in between:\t%d\t%d\n--------------------------------------\n%s\n-------------------------------------------------------\n",max(pt_start - LENFLAG,0), pt_end + LENFLAG, ref_part.c_str());
			al.align(ref_part, contig.data);
			al.dump(fo_full);
			if(al.extract_calls(contig.data, ref_part, reports, fo_vcf, fo_full, fo_vcf_del, fo_full_del, pt.get_reference(), contigSupport, pt_start, pt_end, contigNum, LENFLAG)==0)
			{
				string rc_contig = reverse_complement(contig.data);	
				al.align(ref_part, rc_contig);
				al.dump(fo_full);
				al.extract_calls(rc_contig, ref_part, reports, fo_vcf, fo_full, fo_vcf_del, fo_full_del, pt.get_reference(), contigSupport, pt_start, pt_end, contigNum, LENFLAG);
			}
			contigNum++;
		}
		print_calls(reports, fo_vcf, fo_full, pt.get_cluster_id());
		if((reports.size() == 0 || reports.size()>1) && hybrid == 1)
		{
			reports.clear();
			string outofsga = assemble_with_sga(out_vcf, p, read_length);

			FILE *fcontig = fopen(outofsga.c_str(),"r");
			vector<string> contigs;
			char *line = new char[1000000];
			int contigNum=0;
			int contigSupport =p.size();
			while(fgets(line,1000000,fcontig)!=NULL)
			{
				fgets(line,1000000,fcontig);
				line[strlen(line)-1]='\0';
				string contig = string(line);
				int con_len = contig.length();
				if(check_AT_GC(contig, MAX_AT_GC)==0 || contigSupport <=1 || con_len > max_len + 400 || (pt_end +1000 - (pt_start - 1000)) > MAX_REF_LEN ) continue;
				contigNum++;
				
				fprintf(fo_full,"#################################################################################################\nNEWCONTIG %d IN PARTITION %d length %lu; readNum %d; contig range %d\t%d\n%s\nSUPPORTIVE READS ARE:\n", contigNum, cluster_id, contig.length(), contigSupport, pt_start, pt_end, contig.c_str());
				string ref_part = ref.extract(pt.get_reference(), pt_start - LENFLAG, pt_end + LENFLAG);
				int matrix_size = ref_part.size();
				if(con_len > matrix_size) matrix_size = con_len;
				aligner al(ANCHOR_SIZE, matrix_size);
				fprintf(fo_full, "REF region is in between:\t%d\t%d\n--------------------------------------\n%s\n-------------------------------------------------------\n",max(pt_start - LENFLAG,0), pt_end + LENFLAG,ref_part.c_str());
				al.align(ref_part, contig);
				al.dump(fo_full);
				if(al.extract_calls(contig, ref_part, reports, fo_vcf, fo_full, fo_vcf_del, fo_full_del, pt.get_reference(), contigSupport, pt_start, pt_end, contigNum, LENFLAG)==0)
				{
					string rc_contig = reverse_complement(contig);	
					al.align(ref_part, rc_contig);
					al.dump(fo_full);
					al.extract_calls(rc_contig, ref_part, reports, fo_vcf, fo_full, fo_vcf_del, fo_full_del, pt.get_reference(), contigSupport, pt_start, pt_end, contigNum, LENFLAG);
				}
			}
			print_calls( reports, fo_vcf, fo_full, pt.get_cluster_id());
			reports.clear();
		}
	}	
	fclose(fo_vcf);
	fclose(fo_full);
	fclose(fo_full_del);
	fclose(fo_vcf_del);
}
/*****************************************************/
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
/*************************************************************/
void wo (FILE *f, char *n, char *s, char *q) {
	int m = min(strlen(q), strlen(s));
	s[m] = q[m] = 0;
	fprintf(f, "@%s\n%s\n+\n%s\n", n, s, q);
}
/***********************************************************/
char checkNs (char *s) {
	int nc = 0, c = 0;
	while (*s) nc += (*s == 'N'), s++, c++;
	return (nc > 10 ? 0 : 1);
}
/*********************************************************************************************/
int main(int argc, char **argv)
{

	try {
		if (argc < 2) throw "Usage:\tsniper [mode=(?)]";

		string mode = argv[1];
		if (mode == "verify_sam") {
			if (argc < 7) throw "Usage:\tsniper fastq [sam-file] [output] outputtype oea? orphan?";
			extractor ext(argv[2], argv[3], atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
		}
		else if (mode == "remove_concordant") {
			if (argc != 7) throw "Usage:\tsniper oea [sam-file] [output] outputtype oea? orphan?";
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
			if (argc != 10) throw "Usage:10 parameters needed\tsniper assemble [partition-file] [reference] [range] [output-file-vcf] [output-file-full] [max-len] [read-length] [hybrid]"; 
			assemble(argv[2], argv[3], argv[4], argv[5], argv[6], atoi(argv[7]), atoi(argv[8]), atoi(argv[9]));
		}
		else if (mode == "assemble_orphan_with_sga") {
			if (argc != 3) throw "Usage:3 parameters needed\tsniper assemble_orphan_with_sga [orphan-fastq-file]"; 
			assemble_orphan_with_sga(argv[2]);
		}
		else if (mode == "get_cluster") {
			if (argc != 4) throw "Usage:\tsniper get_cluster [partition-file] [range]";
			genome_partition pt;
			pt.output_partition( argv[2], argv[3]);
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
