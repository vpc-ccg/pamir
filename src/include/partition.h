#ifndef __PARTITION__
#define __PARTITION__

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <zlib.h>
#include <unordered_set>
using namespace std;

class genome_partition {
private:
    FILE* partition_out_file;
    FILE* partition_out_index_file;
    FILE* partition_out_count_file;

    int partition_count = 0;
    int partition_id = 0;

	int distance;

//    int fc;
    int p_start;
    int p_end;
    string p_ref;

    //TODO: more todo
    int start;
    int end;
    FILE* partition_file;





	unordered_map<string, string> oea_mate;
	map<string, vector<vector<string> > > oea2contig;
	map<string, string> map_cont;
	// insert orphan into clusters

	// Orphan contig stats specific to current cluster
	map<string, vector<int> > orphan_contig_stats;
    //TODO: properly rename it to current_cluster;
    vector<pair<pair<string, string>, pair<int,int> > > current_cluster;


	gzFile fp;
	char prev_string[2000];

	const int INSSIZE=450;

private:
	void add_read (string, int, int);
    void update_clusters_with_orphan_contigs ( );
    void dump();
    bool has_next (void);
public:
	genome_partition (const string&, const string&);
	genome_partition (int, const string&,  const string&, const string&, const string&, const string &);
	~genome_partition (void);
	genome_partition(const string&);



    void cluster_reads ();
	int load_orphan( const string &orphan_contig, const string &oea2orphan);
	int load_oea_mates (const string &mate_file) ;

	void get_next (void);
	vector<pair<pair<string, string>, pair<int,int> > > read_partition ();




	int get_start (void);
	int get_end (void);
	string get_reference (void);
	void output_partitions();
	int get_id ();
};


#endif
