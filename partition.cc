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
}

genome_partition::~genome_partition ()
{
	if (fp) gzclose(fp);
}

bool genome_partition::add_read (string read_name, int flag, int loc)
{
	if (read_name.size() > 2 && read_name[read_name.size() - 2] == '/')
		read_name = read_name.substr(0, read_name.size() - 2);
	
	if (read_cache.find(read_name) != read_cache.end()) return false;
	auto it = oea_mate.find(read_name);
	if (it != oea_mate.end()) {
		if (flag & 0x10)
			comp.push_back({{read_name + "+", it->second}, {1, loc}});
		else
			comp.push_back({{read_name + "-", reverse_complement(it->second)},{1, loc}});
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
	comp.clear();
	read_cache.clear();

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
	
	return comp;
}

size_t genome_partition::dump (const vector<pair<pair<string, string>, pair<int,int>>> &vec, FILE *fo, int fc)
{
	size_t pos = ftell(fo);
	fprintf(fo, "%d %d %d %d %s\n", fc, vec.size(), p_start, p_end, p_ref.c_str());
	for (auto &i: vec)
		fprintf(fo, "%s %s %d %d\n", i.first.first.c_str(), i.first.second.c_str(), i.second.first, i.second.second);
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
	const int MAXB = 259072;
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
// To get cluster id t, please specify t-t+1, 
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
			end = tok ? atol(tok) : -1;
		}
		free(dup);
		fprintf(stdout, "extraction [%ld, %ld]\n", start, end-1);

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
	const int MAXB = 8096;
	char pref[MAXB];
	char name[MAXB], read[MAXB];
	string c_file = range + ".cluster";
	fo = fopen(c_file.c_str(), "w");
	foidx = fopen((c_file+".idx").c_str(),"wb");
	fclose(fo);
	fclose(foidx);
	size_t foidx_size;
reset:
//	assert(start < offsets.size());
	if (start >= offsets.size() || start >= end)
		return 0;
	//fprintf(stderr,"Seeking to %d--%d (%lu)\n", start, end, offsets[start]);

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
		//result.push_back({{string(name), string(read)}, loc});
	}
	fclose(fi);
	fclose(fo);
	if ( num_read == 0)
		goto reset;
	fclose(foidx);
	return num_cluster;
}

