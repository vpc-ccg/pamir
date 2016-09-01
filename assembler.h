#ifndef __ASSEMBLER__
#define __ASSEMBLER__

#include<iostream>
#include<string>
#include<vector>
#include<utility>
#include<set>

using namespace std;

struct contig {
	string data;

	struct contig_read {
		int location;
		int in_genome_location;
		string name;
		string data;
	};
	vector<contig_read> read_information;

	contig (void) {
	}
	int supportVal=0;
	int support (void) const { 
		return supportVal;
		//return read_information.size();
	}
	int readInfoSize (void) const { 
		return read_information.size();
	}

	int start = -1;
	int get_start () {
		if (start > -1) return start;
		start=read_information[0].in_genome_location;		
		for (int i=0; i < read_information.size(); i++)
		{
			//printf("ALL:\t%d\t%d\t%d\n",i,start,read_information[i].in_genome_location);
			if(read_information[i].in_genome_location!=-1)
			{
				if(start==-1)
				{
					start=read_information[i].in_genome_location;
				}
				if(start>read_information[i].in_genome_location)
				{
					start = read_information[i].in_genome_location;
				}
			}
		}
		return start;
	}

	int end = -1;
	int get_end () {
		if (end > -1) return end;
		end = read_information[0].in_genome_location + read_information[0].data.size();
		for (int i = 1; i < read_information.size(); i++)
			end = max(end, read_information[i].in_genome_location + (int)read_information[i].data.size());
		return end;
	}
};

class assembler {
private:
	struct sread {
		// first one real loc, second one in-contig loc
		vector<pair<pair<pair<int, int>, pair<string, string>>,int>> reads;
		string data;
		int used;
		sread (void): 
			used(0), data("") {}
	};

private:
	const int max_contig_size;
	const int min_glue_size;
//	static const int prime_seed=157;
	static const int prime_seed=191;

	set<sread*> **hash_p;
	set<sread*> **hash_s;
	set<sread*> reads;

private:
	void initialize (int mcs, int ps);
	void update (sread *r, char mode);
	bool full_compare (const string &a, const string &b, int glue_sz);
	bool assemble_single (int sz, int ha);
	void assemble (void);
	void print ();

public:
	assembler (int, int);
	~assembler (void);

public:
	vector<contig> assemble (const vector<pair<pair<string, string>, pair<int, int>>> &input);
};

#endif
