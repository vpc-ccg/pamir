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

genome_partition::genome_partition (const string & out_prefix){
    partition_out_file = fopen(out_prefix.c_str(), "wb");
    if (partition_out_file == NULL)
        throw("[Genome Partition] Cannot open the file.");

    partition_out_index_file = fopen((out_prefix + ".idx").c_str(), "wb");
    if (partition_out_index_file == NULL)
        throw("[Genome Partition] Cannot open the file.");

    partition_out_count_file = fopen((out_prefix + ".count").c_str(), "w");
    if (partition_out_count_file == NULL)
        throw("[Genome Partition] Cannot open the file.");
}

genome_partition::genome_partition (const string &partition_file_path, const string &range):genome_partition(range) {
    // extracting range [start,end]
    size_t pos;
    start = stoi (range, &pos);
    if (pos < range.size()) {
        end = stoi (range.substr(pos+1));
    } else {
        end = start;
    }

    FILE *partition_file_index = fopen((partition_file_path + ".idx").c_str(), "rb");
    if ( partition_file_index == NULL)
        throw("[Genome Partition] Cannot open the file.");

    Logger::instance().info ("Loading the index file.\n");
    vector<size_t> offsets;
    size_t offset;
    while (fread(&offset, 1, sizeof(size_t), partition_file_index) == sizeof(size_t))
    {
        offsets.push_back(offset);
    }

    fclose(partition_file_index);

    if (start < 1)
        start = 1;
    if ( end > offsets.size()+1 )
        end = offsets.size()+1;

    partition_file = fopen(partition_file_path.c_str(), "rb");
    fseek(partition_file, offsets[start-1], SEEK_SET);
}

genome_partition::genome_partition(int dist, const string &filename, const string &contig_file, const string &oea2orphan, const string& mate_file ,const string &output):genome_partition(output) {
    load_orphan( contig_file, oea2orphan);
    load_oea_mates (mate_file);


    distance = dist;

    fp = gzopen(filename.c_str(), "r");
    gzgets(fp, prev_string, 2000);
    while (!strncmp("@", prev_string, 1)) {
        gzgets(fp, prev_string, 2000);
    }
}

void genome_partition::partify() {
    fc=1;
    while (has_next()) {
        auto p = get_next();
        cluster_count++;
        size_t i = dump(p, partition_out_file, fc);
//        fwrite(&i, 1, sizeof(size_t), partition_out_index_file);
        new_dump();
        fc++;
    }
}


genome_partition::~genome_partition ()
{
    fprintf(partition_out_count_file, "%d\n", cluster_count);

	if (fp) gzclose(fp);

	if (partition_out_file != NULL)
	    fclose (partition_out_file);

    if (partition_out_index_file != NULL)
        fclose (partition_out_index_file);

    if (partition_out_count_file != NULL)
        fclose (partition_out_count_file);

}

int genome_partition::load_oea_mates ( const string &mate_file) {
    FILE *fin = fopen(mate_file.c_str(), "r");
    char name[MAXB], read[MAXB], tmp[MAXB];
    while (fgets(name, MAXB, fin)) {
        fgets(read, MAXB, fin);
        fgets(tmp, MAXB, fin);
        fgets(tmp, MAXB, fin);
        if (strlen(name)>2 && name[strlen(name) - 3] == '/')
            name[strlen(name)-3]='\0';
        read[strlen(read)-1]='\0';
        oea_mate[string(name+1)] = read;
    }
    fclose(fin);
    return 0;
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

void genome_partition::new_dump() {
    size_t pos = ftell(partition_out_file);
    fwrite(&pos, 1, sizeof(size_t), partition_out_index_file);
    fprintf(partition_out_file, "%d %lu %d %d %s\n", fc, current_cluster.size(), p_start, p_end, p_ref.c_str());
    for (auto &i: current_cluster)
        fprintf(partition_out_file, "%s %s %d %d\n", i.first.first.c_str(), i.first.second.c_str(), i.second.first, i.second.second);
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
//	size_t pos = ftell(fo);
//	fprintf(fo, "%d %d %d %d %s\n", fc, nReads + acceptedContigNum, p_start, p_end, p_ref.c_str());
//	for (auto &i: vec)
//		fprintf(fo, "%s %s %d %d\n", i.first.first.c_str(), i.first.second.c_str(), i.second.first, i.second.second);
	if ( acceptedContigNum )
		fprintf(fo, "%s", orphan_info.c_str() );
	size_t pos = 0;
	return pos;
}
int genome_partition::get_cluster_id ()
{
	return fc;
}

//read next partition
vector<pair<pair<string, string>, pair<int,int>>> genome_partition::read_partition ()
{

	int sz, i;
//	vector<pair<pair<string, string>, pair<int,int>>> result;

	char pref[MAXB];
	char name[MAXB], read[MAXB];

//    current_cluster.clear();

reset:
	if (start > end)
		return vector<pair<pair<string, string>, pair<int,int>>>();
    cluster_count++;
	start++;
	fscanf(partition_file, "%d %d %d %d %s\n", &fc, &sz, &p_start, &p_end, pref);
	p_ref = pref;
	current_cluster.resize(0);
    current_cluster.reserve(sz);

	for (i = 0; i < sz; i++) {
		fgets(pref, MAXB, partition_file);
		int loc;
		int support;
		sscanf(pref, "%s %s %d %d", name, read, &support, &loc);
		string sname(name);
		string sread(read);
		current_cluster.push_back({{sname, sread}, {support, loc}});
	}

	if (current_cluster.size() == 0)
		goto reset;
	return current_cluster;
}

void genome_partition::output_partitions() {
    vector < pair < pair < string, string >, pair < int, int >> > p;
    while (1) {
        p = read_partition();
        if (p.size() == 0)
            break;
        new_dump();

//        Logger::instance().info("%d %d %d %d %s\n", get_cluster_id(), p.size(), p_start, p_end, p_ref.c_str());
//        for (auto it = p.begin(); it != p.end(); it++) {
//            Logger::instance().info("%s %s %d %d\n", it->first.first.c_str(), it->first.second.c_str(),
//                                    it->second.first, it->second.second);
//        }

    }
}

