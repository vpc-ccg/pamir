#ifndef __PARTITION_HYBRID__
#define __PARTITION_HYBRID__

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <zlib.h>
#include <unordered_set>
using namespace std;

class genome_partition_hybrid {
private:
    FILE* partition_out_file;
    FILE* partition_out_index_file;
    FILE* partition_out_count_file;

    string map_file;

    int partition_count = 0;
    int partition_id = 0;

	int distance;

    int p_start;
    int p_end;
    string p_ref;

    //TODO: more todo
    int start;
    int end;
    FILE* partition_file;

    int total;

    vector<pair<pair<string, string>, pair<int,int> > > current_cluster;

	char prev_string[2000];

	const int INSSIZE=450;

private:
	void add_read(string, int, int);

public:
	genome_partition_hybrid(const string&, const string&, bool write_to_file = false);
	genome_partition_hybrid(int, const string&,  const string&, const string&, const string&, const string &);
	~genome_partition_hybrid(void);
	genome_partition_hybrid(const string&, bool);

    void cluster_reads(string map_path, int len);

	vector<pair<pair<string, string>, pair<int,int> > > read_partition();
	
	int get_start(void);
	int get_end(void);
	string get_reference(void);
	void output_partitions();
	int get_id();
	int get_total();
};


#endif
