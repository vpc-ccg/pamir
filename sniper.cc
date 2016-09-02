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
//#include "IntervalTree.h"
#include <DeeZ.h>
#include <Sort.h>
//#include "../dz/DeeZ.h"
//#include "../dz/Sort.h"
//#include "../dz/Parsers/SAMParser.h"
//#include "../dz/Parsers/BAMParser.h"
using namespace std;


inline string space (int i) 
{
	return string(i, ' ');
}

inline string itoa (int i)
{
	char c[50];
	sprintf(c, "%d", i);
	return string(c);
}
/******************************************************************/
int check_AT_GC(const string &contig, const double &MAX_AT_GC)
{
	double AT_count = 0, GC_count = 0;
	int clen        = contig.length();
	for(int i = 0; i < clen-1; i++)
	{
		if((contig[i] == 'A' && contig[i+1] == 'T') ||(contig[i] == 'T' && contig[i+1] == 'A')) AT_count++;
		if((contig[i] == 'C' && contig[i+1] == 'G') ||(contig[i] == 'G' && contig[i+1] == 'C'))	GC_count++;
	}
	if(AT_count/(double)clen >= MAX_AT_GC || GC_count/(double)clen >= MAX_AT_GC) 
		return 0;
	else
		return 1;
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
int map_short_contig_to_ref(const string &contig, const string &ref_part, aligner &al,  vector<tuple<string, int, int, string, int>> &reports, FILE *fo_vcf, FILE *fo_full, FILE *fo_vcf_del, FILE *fo_full_del, genome_partition &pt, const int &read_length, const int &contigSupport, const int &contig_start, const int &contig_end, const int &contigNum)
{
	int mapped=0;
	string insertion_content;
	double fwdIden = 0;
	double secFwdIden = 0;
	int fwdS, fwdE, fwdAS, fwdAE, secFwdS, secFwdE, secFwdAS, secFwdAE;
	int insertion_start_loc=-1;
	int insertion_end_loc=-1;
	int emptyInsertion=0;

	int LENFLAG=1000;
//	int LENFLAG=4*read_length;
	al.align(ref_part, contig);
	al.dump(fo_full);
	fwdIden = al.get_identity();
	if(fwdIden>1.0 && (al.get_left_anchor()>5 && al.get_right_anchor()>5) && ( al.get_left_anchor() > 16 || al.get_right_anchor()>16) )
	{
		fprintf(fo_full,"SHORT CONTIG\n");
		if(contig_start-LENFLAG<=0)
			fwdS=al.get_start();
		else
			fwdS = contig_start - LENFLAG + al.get_start();
		if(contig_start-LENFLAG<=0)
			fwdE=al.get_end();
		else
			fwdE = contig_start - LENFLAG + al.get_end();
		fwdAS = 0;
		fwdAE = contig.length()-1;
		fprintf(fo_full,"fwdIden: %f\nfwdAE: %d\tfwdAS: %d\tfwdE: %d\tfwdS: %d\tal.get_start(): %d\n",fwdIden, fwdAE,fwdAS,fwdE,fwdS,al.get_start());			
	
		string rpart=al.get_a();
		string cpart=al.get_b();
		fprintf(fo_full,"REF PART OF MAPPING: %s\n",rpart.c_str());
		int del_len=0;
		int ins_len=0;
		int deletion_start_loc=0;
		int deletion_end_loc=0;;
		string deletion_content;
		int isdeletion=0;
		if(cpart.length()>0)
		{
			int p=0;
			while(p<cpart.length())
			{
				del_len=0;
				for(p;p<cpart.length();p++)
				{
					if(cpart[p]=='-')
					{
						deletion_start_loc=p;
						break;
						ins_len=1;
					}
				}
				for(p;p<cpart.length();p++)
				{
					if(cpart[p]!='-')
					{
						deletion_end_loc=p;
						break;
					}
					else
						del_len++;
				}
				if(del_len>=5&&deletion_start_loc!=-1&&deletion_start_loc+1!=insertion_end_loc)
				{
					isdeletion=1;
					fprintf(fo_full_del,"CONTIG ITSELF:%s\n",contig.c_str());
					fprintf(fo_full_del,"del_start:%d\tdel_end:%d\tdel_end-del_start:%d\n",deletion_start_loc,deletion_end_loc,deletion_end_loc-deletion_start_loc);
					string mapped_ref=al.get_a();
					deletion_content=mapped_ref.substr(deletion_start_loc,deletion_end_loc-deletion_start_loc);
					deletion_start_loc=deletion_start_loc+fwdS;
					if(deletion_content.length()>0)
					{
						fprintf(fo_vcf_del,"%s\t%d\t%d\t%s\t%d\n",pt.get_reference().c_str(),deletion_start_loc,deletion_content.length(),deletion_content.c_str(),contigSupport);											
					}
				}
			}
		}
		if(isdeletion==0 and rpart.length()>0)
		{
			int p=0;
			int decreaseFromInsStartLoc=0;
			while(p<rpart.length())
			{
				ins_len=0;
				for(p;p<rpart.length();p++)
				{
					if(rpart[p]=='-')
					{
						insertion_start_loc=p;
						break;
						ins_len=1;
					}
				}
				for(p;p<rpart.length();p++)
				{
					if(rpart[p]!='-')
					{
						insertion_end_loc=p;
						break;
					}
					else
						ins_len++;
				}
				if(ins_len>=5&&insertion_start_loc!=-1&&insertion_start_loc+1!=insertion_end_loc)
				{
					fprintf(fo_full,"CONTIG ITSELF:%s\n",contig.c_str());
					fprintf(fo_full,"ins_start:%d\tins_end:%d\tins_end-ins_start:%d\n",insertion_start_loc,insertion_end_loc,insertion_end_loc-insertion_start_loc);
					string mapped_con=al.get_b();
					insertion_content=mapped_con.substr(insertion_start_loc,insertion_end_loc-insertion_start_loc);
					insertion_start_loc=insertion_start_loc+fwdS;
					if(insertion_content.length()>0)
					{
						fprintf(fo_full,"%s\t%d\t%d\t%s\n------------------------------------------------------\n",pt.get_reference().c_str(),insertion_start_loc,insertion_content.length(),insertion_content.c_str());
						//fprintf(fo_vcf,"%s\t%d\t%d\t%s\t%d\n",pt.get_reference().c_str(),insertion_start_loc,insertion_content.length(),insertion_content.c_str(),contigSupport);					
						reports.push_back(tuple<string, int, int, string, int>(pt.get_reference(), insertion_start_loc, insertion_content.length(), insertion_content,contigSupport ) );
						mapped=1;
					}
				}
				else
				{
					decreaseFromInsStartLoc+=ins_len;
					fprintf(fo_full,"INSERTION LENGTH IS NOT LONG ENOUGH!\n------------------------------------------------------\n");
				}
				p++;
			}
		}

	}
	else
	{
		fprintf(fo_full,"fwdIden: %f is less than 20 so no insertion from this contig.\n",fwdIden);
	}
	return mapped;
}
/********************************************************************************/
void evaluate_contig(const string &contig, vector<tuple<string, int, int, string, int>> &reports,  genome &ref, genome_partition &pt, FILE *fo_vcf, FILE *fo_full, FILE *fo_vcf_del, FILE *fo_full_del, const int &read_length, const int &pt_start, const int &pt_end, const int &contigNum, const int &contigSupport, const int &pCount, const int &ANCHOR_SIZE )
{
	fprintf(fo_full,"#################################################################################################\nNEWCONTIG %d IN PARTITION %d length %d; readNum %d; contig range %d\t%d\n%s\n",contigNum,pCount,contig.length(),contigSupport,pt_start,pt_end,contig.c_str());
	
	//if contig length is greater than 2 x read_length than map its prefix and suffix onto the reference section
//	if(contig.length()>200000)
/*	if(contig.length()>2*read_length)
	{
		string ref_part = ref.extract(pt.get_reference(), pt_start-1000, pt_end+1000);
		fprintf(fo_full, "contig_start: %d\tcontig_end: %d",pt_start-1000, pt_end+1000);
		aligner al(ANCHOR_SIZE,ref_part.length());
		fprintf(fo_full, "REF region is in between:\t%d\t%d\n--------------------------------------\n%s\n-------------------------------------------------------\n",max(pt_start- 1000,0), pt_end + 1000,ref_part.c_str());
		if(map_to_ref(contig, ref_part, al, reports, fo_vcf, fo_full, fo_vcf_del, fo_full_del, pt, read_length, contigSupport, pt_start, pt_end, contigNum)==0)
		{	
			char *rc_contig=new char[contig.length()+1];
			char *c_contig=new char[contig.length()+1];
			strcpy(c_contig, contig.c_str());
			copy_string_rc(c_contig,rc_contig);
			string src_contig=string(rc_contig);
			map_to_ref(src_contig, ref_part, al, reports, fo_vcf, fo_full, fo_vcf_del, fo_full_del, pt, read_length, contigSupport, pt_start, pt_end, contigNum);
		}
	}
*/	//if(!(con_len> 2*read_length)) map the contig itself onto the reference section.
//	else
	{
//		int LENFLAG         = 4*read_length;
		int LENFLAG         = 1000;
		int start_point = pt_start-LENFLAG;
		int end_point   = pt_end+LENFLAG;
		string ref_part = ref.extract(pt.get_reference(), start_point, end_point);
		fprintf(fo_full, "contig_start: %d\tcontig_end: %d\n",start_point, end_point);
		int matrix_size = ref_part.size();
		if(contig.length() > ref_part.size())
			matrix_size = contig.length();
		aligner al(ANCHOR_SIZE,matrix_size);
		fprintf(fo_full, "REF region is in between:\t%d\t%d\n--------------------------------------\n%s\n-------------------------------------------------------\n",max(pt_start - LENFLAG,0), pt_end + LENFLAG,ref_part.c_str());
		if(map_short_contig_to_ref(contig, ref_part, al, reports, fo_vcf, fo_full, fo_vcf_del, fo_full_del, pt, read_length, contigSupport, pt_start, pt_end, contigNum)==0)
		{
			string rc_contig = reverse_complement(contig);	
			map_short_contig_to_ref(rc_contig,ref_part, al, reports, fo_vcf, fo_full, fo_vcf_del, fo_full_del, pt, read_length, contigSupport, pt_start, pt_end, contigNum);
		}
	}
}
/****************************************************************/
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
/**************************************************************************************/
void assemble (const string &partition_file, const string &reference, const string &range, const string &out_vcf, const string &out_full, int max_len, int read_length, const int &hybrid)
{
	const double MAX_AT_GC = 0.7;
	const int ANCHOR_SIZE = 16;
	const int MAX_REF_LEN=300000000;
	assembler as(max_len, 15);
	genome ref(reference.c_str());
	genome_partition pt;
	string con_pre,con_suf;
	int coverage[(max_len + 200)];
	double fwdIden = 0;
	double secFwdIden = 0;
	int fwdS, fwdE, fwdAS, fwdAE, secFwdS, secFwdE, secFwdAS, secFwdAE;
	FILE *fo_vcf = fopen(out_vcf.c_str(), "w");
	FILE *fo_vcf_del = fopen((out_vcf+"_del").c_str(), "w");
	FILE *fo_full = fopen(out_full.c_str(), "w");	
	FILE *fo_full_del = fopen((out_full+".del").c_str(), "w");	
	
	string insertion_content;
	int insertion_start_loc=-1;
	int insertion_end_loc=-1;
	int pCount=0;
	int LENFLAG=250;
	int contig_start=-1;
	int contig_end=-1;
	int pt_start=-1;
	int pt_end=-1;
	while (1) 
	{
		auto p = pt.read_partition(partition_file, range);
		pCount=pCount+1;
		pt_start=pt.get_start();
		pt_end=pt.get_end();
		if (!p.size()) 
			break;
		if (p.size() > 100000) 	continue;
		if(p.size()<=2)	continue;	
		auto contigs = as.assemble(p);
		fprintf(fo_full,"PARTITION %d | read number in partition= %d; contig number in partition= %d; partition range=%d\t%d\n*******************************************************************************************\n",pCount,p.size(),contigs.size(),pt_start, pt_end);
					
		int contigNum=1;
		int ATGCRichContigNum=0;
		vector< tuple< string, int, int, string, int > > reports;
		for (auto &contig: contigs)
		{
			contig_start=contig.get_start();
			contig_end=contig.get_end();
			int contigSupport=contig.support();
			if(contig_start<0)
			{
				contig_start=pt_start;
				contig_end=pt_end;
			}
			int con_len = contig.data.length();
			if(!check_AT_GC(contig.data, MAX_AT_GC)) continue;
			
			if(contigSupport > 1)
			{				
				if (con_len > max_len + 400 || contig_end+1000-(contig_start-1000) > MAX_REF_LEN) 
				{
					fprintf(stderr, "WEIRD CASE Too long region -- Contig %d Genomic %d!\n", contig.data.size(), contig_end+1000-(contig_start-1000));

					for (auto &e: contig.read_information)
					{
						fprintf(stderr, "%d %d %s %s\n", e.in_genome_location, e.location, e.data.c_str(), e.name.c_str());
					}
					continue;
				}			
				fprintf(fo_full,"#################################################################################################\nNEWCONTIG %d IN PARTITION %d length %d; readNum %d; contig range %d\t%d\n%s\nSUPPORTIVE READS ARE:\n",contigNum,pCount,contig.data.size(),contigSupport,contig_start,contig_end,contig.data.c_str());
				for(int z=0;z<contig.read_information.size();z++)
					fprintf(fo_full,"%s %d %d %s\n",contig.read_information[z].name.c_str(),contig.read_information[z].location,contig.read_information[z].in_genome_location, contig.read_information[z].data.c_str());		

	/*			for (int k=0; k < con_len; k++) coverage[k]=0;
				for (int j = 0; j < contig.readInfoSize(); j++)
				{
					for (int k=0; k < contig.read_information[j].data.length(); k++)
					{
						if(k+contig.read_information[j].location >= -1)
						{
							coverage[k+contig.read_information[j].location]++;
	
						}
					}
				}*/
				evaluate_contig(contig.data, reports, ref, pt, fo_vcf, fo_full, fo_vcf_del, fo_full_del, read_length, pt_start, pt_end, contigNum, contigSupport, pCount, ANCHOR_SIZE );
				contigNum++;
			}
		}
		int difflocorlength =0;
		if(reports.size()>0)
		{	
			fprintf(fo_full,"DECIDED BY INHOUSE ASSEMBLER\n");
			for(int r=0;r<reports.size();r++)
			{
				fprintf(fo_full,"%s\t%d\t%d\t%s\t%d\n------------------------------------------------------\n",get<0>(reports[r]).c_str(), get<1>(reports[r]), get<2>(reports[r]), get<3>(reports[r]).c_str(),get<4>(reports[r]));
				fprintf(fo_vcf,"INHOUSE\t%s\t%d\t%d\t%s\t%d\n", get<0>(reports[r]).c_str(), get<1>(reports[r]), get<2>(reports[r]), get<3>(reports[r]).c_str(), get<4>(reports[r]));					
			}
		}
		else if((reports.size() == 0 || reports.size()>1) && hybrid == 1)
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
				contigNum++;
				if(!check_AT_GC(contig,MAX_AT_GC)) continue;
				if (con_len > max_len + 400 || pt_end+1000-(pt_start-1000) > MAX_REF_LEN) 
				{
					fprintf(stderr, "WEIRD CASE Too long region -- Contig %d Genomic %d!\n", contig.length(), pt_end+1000-(pt_start-1000));
					continue;
				}			
				fprintf(fo_full,"DECIDED BY SGA ASSEMBLER\n");
				evaluate_contig(contig, reports, ref, pt, fo_vcf, fo_full, fo_vcf_del, fo_full_del, read_length, pt_start, pt_end, contigNum, contigSupport, pCount, ANCHOR_SIZE );
				for(int r =0;r<reports.size();r++)
					fprintf(fo_vcf,"SGA\t%s\t%d\t%d\t%s\t%d\n", get<0>(reports[r]).c_str(), get<1>(reports[r]), get<2>(reports[r]), get<3>(reports[r]).c_str(), get<4>(reports[r]));					
			}
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


void removeUnmapped (const string &path, const string &result)
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

void sortSAM (const string &path, const string &result) 
{
	sortFile(path, result, 2 * GB);
}

void wo (FILE *f, char *n, char *s, char *q) {
	int m = min(strlen(q), strlen(s));
	s[m] = q[m] = 0;
	fprintf(f, "@%s\n%s\n+\n%s\n", n, s, q);
}

char checkNs (char *s) {
	int nc = 0, c = 0;
	while (*s) nc += (*s == 'N'), s++, c++;
	return (nc > 10 ? 0 : 1);
}

void copy_string_reverse( char *src, char *dest)
{
	int limit=strlen(src);
	for( int i = 0; i < limit; i ++)
	{
		dest[i]=src[limit-1-i];
	}
	dest[limit]='\0';
}

/*********************************************************************************************/
int main(int argc, char **argv)
{

	try {
		if (argc < 2) throw "Usage:\tsniper [mode=(?)]";

		string mode = argv[1];
		if (mode == "verify_sam") {
			if (argc < 4) throw "Usage:\tsniper fastq [sam-file] [output]";
			extractor ext(argv[2], argv[3],3,1,1);
		}
		else if (mode == "remove_concordant") {
			if (argc != 4) throw "Usage:\tsniper oea [sam-file] [output]";
			extractor ext(argv[2], argv[3],2,1,1);
		}
		else if (mode == "mask" || mode == "maski") {
			if (argc != 6) throw "Usage:\tsniper mask/maski [repeat-file] [reference] [output] [padding]";
			mask(argv[2], argv[3], argv[4], atoi(argv[5]), mode == "maski");
		}
		else if (mode == "sort") {
			if (argc != 4) throw "Usage:\tsniper sort [sam-file] [output]";
			sortSAM(argv[2], argv[3]);
		}
		else if (mode == "rm_unmap") {
			if (argc != 6) throw "Usage:\tsniper rm_unmap [fq-file] [output] outputtype oea? orphan?";
			removeUnmapped(argv[2], argv[3]);
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
