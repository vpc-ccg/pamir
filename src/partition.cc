#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cstdio>
#include <cassert>
#include <zlib.h>
#include <cstring>
#include <cstdlib>
#include "common.h"
#include "partition.h"
#define MAXB  809600

using namespace std;

genome_partition::genome_partition () 
{
	distance = 0;
	fp = 0;
}

genome_partition::genome_partition (const string &filename, int dist, const unordered_map<string, string> &om)
{
	distance = dist;
	auto g=om.begin();
	oea_mate = om;

	fp = gzopen(filename.c_str(), "r");	
	gzgets(fp, prev_string, 2000);
	while(  !strncmp("@", prev_string, 1) )
	{
	gzgets(fp, prev_string, 2000);
	}
}

genome_partition::~genome_partition ()
{
	if (fp) gzclose(fp);
}

int genome_partition::load_orphan( const string &orphan_contig, const string &oea2orphan )
{
	//int INSSIZE=450;
	char rname[MAX_CHAR], gname[MAX_CHAR], cont[MAX_CHAR], readline[MAX_CHAR], tmp[MAX_CHAR];
	string rstr, cstr;
	int flag, pos, qual, npos, tlen, z;
	char del[] = " \t\n";
	int i;

	// loading orphan contig
	FILE *fin = fopen( orphan_contig.c_str(), "r");
	while( NULL != fgets( gname, MAX_CHAR, fin ) )
	{
		fgets( cont, MAX_CHAR, fin);
		i = strcspn(gname, del);
		//rstr = string(gname).substr(1,strlen(gname)-2); // excluding leading char and newline
		rstr = string(gname).substr(1,i-1); // excluding leading char and newline
		cstr = string(cont).substr(0, strlen(cont)-1); // excluding leading char and newline
		map_cont[rstr] = cstr;
	}
	fclose(fin);
    Logger::instance().info("Loading orphan contig file\n");
	
	fin = fopen( oea2orphan.c_str(), "r");
	while( NULL != fgets( readline, MAX_CHAR, fin ) )
	{	
		sscanf(readline, "%s %d %s %d %d %s %s %d %d %s %d\n", rname, &flag, gname,  &pos, &qual, tmp, tmp, &npos, &tlen, cont, &z);
		rstr = string( rname );
		if (rstr.size() > 2 && rstr[rstr.size() - 2] == '/')
			rstr = rstr.substr(0, rstr.size() - 2);

		vector<string> tmpv;
		tmpv.push_back(string(gname));
		tmpv.push_back( ( (flag & 0x10) == 0x10 )?"-":"+" );

		int contiglen = map_cont[gname].length();
		int readlen   = strlen(cont);
		if ( pos < INSSIZE)
			tmpv.push_back("l");
		else if(pos > contiglen-INSSIZE-readlen)
			tmpv.push_back("r");
		else
			tmpv.push_back("w");

		if ( map_token.find(rstr) != map_token.end() ) {
			map_token[rstr].push_back(tmpv);
		}
		else
		{
			vector<vector<string> > tmp_vec;
			tmp_vec.push_back(tmpv);
			map_token[rstr] = tmp_vec;
		}
	}
	fclose(fin);
    Logger::instance().info("Updating OEA contigs\n");
	return 0;
}


bool genome_partition::add_read (string read_name, int flag, int loc)
{
	char sign;
	if (read_name.size() > 2 && read_name[read_name.size() - 2] == '/')
		read_name = read_name.substr(0, read_name.size() - 2);
	
	if (read_cache.find(read_name) != read_cache.end()) return false;
	auto it = oea_mate.find(read_name);
	if (it != oea_mate.end()) {
		if (flag & 0x10)
			{	current_cluster.push_back({{read_name + "+", it->second}, {1, loc}}); sign = '+';}
		else
			{	current_cluster.push_back({{read_name + "-", reverse_complement(it->second)},{1, loc}}); sign= '-';}
		// adjust associated orphans for a cluster
		map<string, vector<int> >::iterator mit;
		if (map_token.find(read_name) != map_token.end()) {
			vector<vector<string> >::iterator nit;
			for (nit=map_token[read_name].begin(); nit != map_token[read_name].end(); nit++){
				if (myset.find((*nit)[0]) == myset.end()) {
					vector<int> tmpev;
					tmpev.push_back(0);
					tmpev.push_back(0);
					tmpev.push_back(0);
					tmpev.push_back(0);
					myset[(*nit)[0]] = tmpev;
					//fprintf( stdout, "F_%s_F\n", (*nit)[0].c_str() );
				}

				if (sign != (*nit)[1][0])
					myset[(*nit)[0]][1]++;
				else
					myset[(*nit)[0]][0]++;

				if ((*nit)[2] == "l")
					myset[(*nit)[0]][2]++;
				else
					myset[(*nit)[0]][3]++;
			}
		}

		return true;
	}
	else return false;
}

bool genome_partition::has_next ()
{
	return !gzeof(fp);
}

int genome_partition::get_start ()
{
	return p_start;
}

int genome_partition::get_end ()
{
	return p_end;
}

string genome_partition::get_reference ()
{
	return p_ref;
}

vector<pair<pair<string, string>, pair<int,int>>> genome_partition::get_next ()
{
	current_cluster.clear();
	read_cache.clear();
	myset.clear();

	char read_name[1000], ref_name[1000];
	int flag, loc;
	int ploc;
	sscanf(prev_string,"%s %d %s %d", read_name, &flag, ref_name, &loc);
	p_ref = ref_name;
	p_start = loc;
	p_end = loc;
	ploc = loc;
	add_read(read_name, flag, loc);
	
	while (has_next()&&gzgets(fp, prev_string, 2000)) {
		sscanf(prev_string,"%s %d %s %d", read_name, &flag, ref_name, &loc);
		if (loc - p_start > distance || ref_name != p_ref)
			break;
		
		add_read(read_name, flag, loc);
		ploc = loc;
		p_end = loc;
	}
	
	return current_cluster;
}

size_t genome_partition::dump (const vector<pair<pair<string, string>, pair<int,int>>> &vec, FILE *fo, int fc)
{
	// Inserting possible orphan contig
	int acceptedContigNum =0;
	int nReads = (int)vec.size();
	int tie= 0;
	string orphan_info;
	orphan_info.reserve(20000);
	
	map<string, vector<int> >::iterator mit;
    Logger::instance().info("CLUSTER ID: %d\n", fc);
	for (mit=myset.begin(); mit!=myset.end(); mit++)
	{
		string revcontent;
		//fprintf(flog,"readNum: %d\tcontigname: %s\tleft: %d\tright: %d\tforward: %d\treverse: %d\n", nReads, (*mit).first.c_str(), myset[(*mit).first][2], myset[(*mit).first][3], myset[(*mit).first][0], myset[(*mit).first][1]);
        Logger::instance().info("readNum: %d\tcontigname: %s\tleft: %d\tright: %d\tforward: %d\treverse: %d\n", nReads, (*mit).first.c_str(), myset[(*mit).first][2], myset[(*mit).first][3], myset[(*mit).first][0], myset[(*mit).first][1]);
		
		int contiglen = map_cont[(*mit).first].length(); 
		if(contiglen <= INSSIZE)
		{
			if( 0 == myset[(*mit).first][2] )
			{
			    Logger::instance().info("Short contig and no left or right flank supporters\n");
				continue;
			}
		}
		else 
		{
			if((myset[(*mit).first][2]==0|| myset[(*mit).first][3]==0) && myset[(*mit).first][2] +myset[(*mit).first][3] < 0.3*nReads )
			{
                Logger::instance().info("either no OEAs mapping on left or right flank and not enough left + right flank support < 30%%\n");
				continue;
			}
		}
		if( myset[(*mit).first][0] ==  myset[(*mit).first][1] )
		{
			orphan_info += (*mit).first + " " + map_cont[(*mit).first] + " 0 -1\n";
			string content = map_cont[(*mit).first];
			revcontent = reverse_complement(content);
			orphan_info += (*mit).first + "_rv " + revcontent+ " 0 -1\n";
			tie++;
			acceptedContigNum +=2;
		}
		else
		{
			if(myset[(*mit).first][0] > myset[(*mit).first][1])
			{
				revcontent = map_cont[(*mit).first];
			}
			else if(myset[(*mit).first][1] > myset[(*mit).first][0])
			{
				string content = map_cont[(*mit).first];
				revcontent = reverse_complement(content);
			}
			orphan_info += (*mit).first + " " + revcontent+ " 0 -1\n";
			acceptedContigNum ++;
		}
	}
	if ( tie ) { Logger::instance().info("Forward and reverse were in tie condition : %d times.\n", tie);}

	// Actual Partiton Content 
	size_t pos = ftell(fo);
	fprintf(fo, "%d %d %d %d %s\n", fc, nReads + acceptedContigNum, p_start, p_end, p_ref.c_str());
	for (auto &i: vec)
		fprintf(fo, "%s %s %d %d\n", i.first.first.c_str(), i.first.second.c_str(), i.second.first, i.second.second);
	if ( acceptedContigNum )
		fprintf(fo, "%s", orphan_info.c_str() );
	return pos;
}
int genome_partition::get_cluster_id ()
{
	return fc;
}

vector<pair<pair<string, string>, pair<int,int>>> genome_partition::read_partition (const string &partition_file, const string &range)
{
	static int start = -1, end = -1;
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
		FILE *fidx = fopen((partition_file + ".idx").c_str(), "rb");
		if(fidx==NULL)
		{
			printf(".idx file does not exist!\n");
			exit(-1);
		}
		else
		{
			size_t offset;
			while (fread(&offset, 1, sizeof(size_t), fidx) == sizeof(size_t))
			{
				offsets.push_back(offset);
			}
		}
		fclose(fidx);
	}

	FILE *fi;
	int sz, i;
	vector<pair<pair<string, string>, pair<int,int>>> result;

	char pref[MAXB];
	char name[MAXB], read[MAXB];

reset:
	//	assert(start < offsets.size());
	if (start >= offsets.size() || start >= end)
		return vector<pair<pair<string, string>, pair<int,int>>>();
	//fprintf(stderr,"Seeking to %d--%d (%lu)\n", start, end, offsets[start]);
	fi = fopen(partition_file.c_str(), "rb");
	fseek(fi, offsets[start++], SEEK_SET);
	fscanf(fi, "%d %d %d %d %s\n", &fc, &sz, &p_start, &p_end, pref);
	p_ref = pref;
	result.resize(0);
	result.reserve(sz);

	for (i = 0; i < sz; i++) {
		fgets(pref, MAXB, fi);
		int loc;
		int support;
		sscanf(pref, "%s %s %d %d", name, read, &support, &loc);
		string sname(name);
		string sread(read);
		result.push_back({{sname, sread}, {support, loc}});
	}
	fclose(fi);
	if (result.size() == 0)
		goto reset;
	return result;
}
////////////////////////// Given partition file, range x-y
// Output Cluster from x to (y-1)  to x-y.cluster
// To get cluster id t, please specify t, 
// otherwise report t to the end of partition
int genome_partition::output_partition (const string &partition_file, const string &range)
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

		FILE *fidx = fopen((partition_file + ".idx").c_str(), "rb");
		size_t offset;
		while (fread(&offset, 1, sizeof(size_t), fidx) == sizeof(size_t))
			offsets.push_back(offset);
		fclose(fidx);
	}
	FILE *fi, *fo, *foidx;
	int sz, i;
	int cluster_id;
	int num_cluster = 0, num_read = 0;

	char pref[MAXB];
	char name[MAXB], read[MAXB];
	string c_file = range + ".cluster";
	fo = fopen(c_file.c_str(), "w");
	foidx = fopen((c_file+".idx").c_str(),"wb");
	fclose(fo);
	fclose(foidx);
	size_t foidx_size;
reset:
	if (start >= offsets.size() || start >= end)
		return 0;
	fi = fopen(partition_file.c_str(), "rb");
	fo = fopen(c_file.c_str(), "a");
	fseek(fi, offsets[start++], SEEK_SET);
	foidx = fopen((c_file+".idx").c_str(),"ab");
	fscanf(fi, "%d %d %d %d %s\n", &cluster_id, &sz, &p_start, &p_end, pref);
	foidx_size = ftell(fo);
	fwrite(&foidx_size,1,sizeof(size_t),foidx);
	fprintf(fo, "%d %d %d %d %s\n", cluster_id, sz, p_start, p_end, pref);
	p_ref = pref;

	num_read = 0;
	for (i = 0; i < sz; i++) {
		fgets(pref, MAXB, fi);
		int loc;
		int support;
		sscanf(pref, "%s %s %d %d", name, read, &support, &loc);
		fprintf( fo, "%s %s %d %d\n", name, read, support, loc);
		num_read = 0;
	}
	fclose(fi);
	fclose(fo);
	fclose(foidx);
	if ( num_read == 0)
		goto reset;
	return num_cluster;
}

