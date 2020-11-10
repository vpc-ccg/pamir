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

genome_partition::genome_partition (const string & out_prefix, bool write_to_files): partition_out_file(NULL), partition_out_index_file(NULL), partition_out_count_file(NULL), fp(NULL){
    if(write_to_files){
        partition_out_file = fopen(out_prefix.c_str(), "wb");
        if (partition_out_file == NULL){
            throw("[Genome Partition] Cannot open the file.");
        }
        partition_out_index_file = fopen((out_prefix + ".idx").c_str(), "wb");
        if (partition_out_index_file == NULL){
            throw("[Genome Partition] Cannot open the file.");
        }
        partition_out_count_file = fopen((out_prefix + ".count").c_str(), "w");
        if (partition_out_count_file == NULL){
            throw("[Genome Partition] Cannot open the file.");
        }
    }
}

genome_partition::genome_partition (const string &partition_file_path, const string &range, bool write_to_file):genome_partition(range, write_to_file) {
    // extracting range [start,end]
    size_t pos;
    start = stoi (range, &pos);
    if (pos < range.size()) {
        end = stoi (range.substr(pos+1));
    } else {
        end = start;
    }

    // reading the index file of partitions
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
    if ( end > offsets.size() )
        end = offsets.size();

    partition_file = fopen(partition_file_path.c_str(), "rb");
    fseek(partition_file, offsets[start-1], SEEK_SET);
}

genome_partition::genome_partition(int dist, const string &filename, const string &contig_file, const string &oea2orphan, const string& mate_file ,const string &output):genome_partition(output, true) {

    load_orphan( contig_file, oea2orphan);
    load_oea_mates (mate_file);


    distance = dist;

    fp = gzopen(filename.c_str(), "r");
    gzgets(fp, prev_string, MAXB);
    while (!strncmp("@", prev_string, 1)) {
        gzgets(fp, prev_string, MAXB);
    }
}

void genome_partition::cluster_reads() {
    partition_id=1;
    while (has_next()) {
        get_next();
        partition_count++;
        update_clusters_with_orphan_contigs();
        dump();
        partition_id++;
    }
}


genome_partition::~genome_partition (){
    if( fp != NULL){
        gzclose(fp);
    }
    if (partition_file != NULL){
        fclose(partition_file);
    }
    if (partition_out_file != NULL){
        fclose (partition_out_file);
    }
    if (partition_out_index_file != NULL){
        fclose (partition_out_index_file);
    }
    if (partition_out_count_file != NULL){
        fprintf(partition_out_count_file, "%d\n", partition_count);
        fclose (partition_out_count_file);
    }
}

int genome_partition::load_oea_mates ( const string &mate_file) {
    gzFile fin = gzopen(mate_file.c_str(), "r");
    char name[MAXB], read[MAXB], tmp[MAXB];
    while (gzgets(fin, name, MAXB)) {
        gzgets(fin, read, MAXB);
        gzgets(fin, tmp, MAXB);
        gzgets(fin, tmp, MAXB);
        if (strlen(name)>2 && name[strlen(name) - 3] == '/')
            name[strlen(name)-3]='\0';
        read[strlen(read)-1]='\0';
        oea_mate[string(name+1)] = read;
    }
    gzclose(fin);
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

        //TODO Fix
        if ( pos < INSSIZE)
			tmpv.push_back("l");
		else if(pos > contiglen-INSSIZE-readlen)
			tmpv.push_back("r");
		else
			tmpv.push_back("w");

		if ( oea2contig.find(rstr) != oea2contig.end() ) {
			oea2contig[rstr].push_back(tmpv);
		}
		else
		{
			vector<vector<string> > tmp_vec;
			tmp_vec.push_back(tmpv);
			oea2contig[rstr] = tmp_vec;
		}
	}
	fclose(fin);
    Logger::instance().info("Updating OEA contigs\n");
	return 0;
}


void genome_partition::add_read(string read_name, int flag, int loc) {
    if (read_name.size() > 2 && read_name[read_name.size() - 2] == '/') {
        read_name = read_name.substr(0, read_name.size() - 2);
    }
    auto it = oea_mate.find(read_name);
    if (it != oea_mate.end()) {
        char sign;
        if (flag & 0x10) {
            // push on the reverse of its mate
            current_cluster.push_back({{read_name + "+", it->second},
                                       {1,               loc}});
            sign = '+';
        } else {
            current_cluster.push_back({{read_name + "-", reverse_complement(it->second)},
                                       {1,               loc}});
            sign = '-';
        }

        auto vec_ptr_of_read = oea2contig.find(read_name);
        if (vec_ptr_of_read != oea2contig.end()) {
            for (auto nit = vec_ptr_of_read->second.begin(); nit != vec_ptr_of_read->second.end(); ++nit) {
                string &ctg_name = (*nit)[0];
                if (orphan_contig_stats.find(ctg_name) == orphan_contig_stats.end()) {
                    orphan_contig_stats[ctg_name] = {0, 0, 0, 0};
                }

                // updating the reverse/forward of mapping per orphan contig
                if (sign == (*nit)[1][0]) {
                    orphan_contig_stats[ctg_name][0]++;
                } else {
                    orphan_contig_stats[ctg_name][1]++;
                }

                // Updating boundry statistic per orphan contig
                //TODO Fix Fix Fix Fix
                if ((*nit)[2] == "l") {
                    orphan_contig_stats[ctg_name][2]++;
                } else {
                    orphan_contig_stats[ctg_name][3]++;
                }
            }
        }
    }
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

void genome_partition::get_next () {
	current_cluster.clear();
	orphan_contig_stats.clear();

	char read_name[MAXB], ref_name[MAXB];
	int flag, cur_loc;
	sscanf(prev_string,"%s %d %s %d", read_name, &flag, ref_name, &cur_loc);
	p_ref = ref_name;
	p_start = cur_loc;
	p_end = cur_loc;
	add_read(read_name, flag, cur_loc);

	while (has_next()&&gzgets(fp, prev_string, MAXB)) {
		sscanf(prev_string,"%s %d %s %d", read_name, &flag, ref_name, &cur_loc);
		if (cur_loc - p_start > distance || ref_name != p_ref)
			break;

		add_read(read_name, flag, cur_loc);
		p_end = cur_loc;
	}
}

void genome_partition::dump() {
	assert (current_cluster.size () > 0);
    size_t pos = ftell(partition_out_file);
    fwrite(&pos, 1, sizeof(size_t), partition_out_index_file);
    fprintf(partition_out_file, "%d %lu %d %d %s\n", partition_id, current_cluster.size(), p_start, p_end, p_ref.c_str());
    for (auto &i: current_cluster) {
        fprintf(partition_out_file, "%s %s %d %d\n", i.first.first.c_str(), i.first.second.c_str(), i.second.first,
                i.second.second);
    }
}

void genome_partition::update_clusters_with_orphan_contigs() {

    size_t tie = 0;
    size_t nReads = current_cluster.size();

    Logger::instance().info("CLUSTER ID: %d\n", partition_id);
    for (auto mit = orphan_contig_stats.begin(); mit != orphan_contig_stats.end(); mit++) {

        Logger::instance().info("readNum: %lu\tcontigname: %s\tleft: %d\tright: %d\tforward: %d\treverse: %d\n", nReads,
                                (*mit).first.c_str(), orphan_contig_stats[(*mit).first][2],
                                orphan_contig_stats[(*mit).first][3], orphan_contig_stats[(*mit).first][0],
                                orphan_contig_stats[(*mit).first][1]);

        int contiglen = map_cont[(*mit).first].length();
        //TODO Fix
        if (contiglen <= INSSIZE) {
            if (0 == orphan_contig_stats[(*mit).first][2]) {
                Logger::instance().info("Short contig and no left flank supporters\n");
                continue;
            }
        }
        else{
            if ((orphan_contig_stats[(*mit).first][2] == 0 || orphan_contig_stats[(*mit).first][3] == 0) &&
                orphan_contig_stats[(*mit).first][2] + orphan_contig_stats[(*mit).first][3] < 0.3 * nReads) {
                Logger::instance().info(
                    "either no OEAs mapping on left or right flank and not enough left + right flank support < 30%%\n");
                continue;
            }
        }
        string content = map_cont[(*mit).first];
        if (orphan_contig_stats[(*mit).first][0] == orphan_contig_stats[(*mit).first][1]) {
            current_cluster.push_back({{mit->first, map_cont[(*mit).first]},
                                       {0,          -1}});
            current_cluster.push_back({{mit->first, reverse_complement(content)},
                                       {0,          -1}});
            ++tie;
        } else {
            if (orphan_contig_stats[(*mit).first][0] > orphan_contig_stats[(*mit).first][1]) {
                current_cluster.push_back({{mit->first, content},
                                           {0,          -1}});
            } else if (orphan_contig_stats[(*mit).first][1] > orphan_contig_stats[(*mit).first][0]) {
                current_cluster.push_back({{mit->first, reverse_complement(content)},
                                           {0,          -1}});
            }
        }
    }
    if (tie > 0) {
        Logger::instance().info("Forward and reverse were in tie condition : %d times.\n", tie);
    }
}

int genome_partition::get_id () {
	return partition_id;
}

//read next partition
vector <pair<pair < string, string>, pair<int, int>>> genome_partition::read_partition() {

	int sz, i;
	char pref[MAXB];
	char name[MAXB], read[MAXB];
//reset:
	if (start > end)
		return vector < pair < pair < string, string >, pair < int, int >> > ();
	partition_count++;
	start++;
	fscanf(partition_file, "%d %d %d %d %s\n", &partition_id, &sz, &p_start, &p_end, pref);
	p_ref = pref;
	current_cluster.resize(0);
	current_cluster.reserve(sz);

	for (i = 0; i < sz; i++) {
		fgets(pref, MAXB, partition_file);
		int loc, support;
		sscanf(pref, "%s %s %d %d", name, read, &support, &loc);
		current_cluster.push_back({{string(name), string(read)},
								   {support, loc}});
	}
//	if (current_cluster.size() == 0)
//		goto reset;
	return current_cluster;
}

void genome_partition::output_partitions() {
    while (1) {
        auto p = read_partition();
        if (p.size() == 0) {
            break;
        }
        dump();
    }
}

