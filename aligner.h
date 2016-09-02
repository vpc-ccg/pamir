#ifndef __ALIGNER__
#define __ALIGNER__

#include<string>

using namespace std;

class aligner
{
private:
	int REF_LEN = 10000;
	int CON_LEN =10000;
	int MAX_SIDE = 15000;
/*	static const int MATCH_SCORE = 20;
	static const int MISMATCH_SCORE = -20;
	static const int GAP_OPENING_SCORE = -100;
	static const int GAP_EXTENSION_SCORE = 0;
	*/
	static const int MATCH_SCORE = 20;
	static const int MISMATCH_SCORE = -1000;
	static const int GAP_OPENING_SCORE = -1000;
	static const int GAP_EXTENSION_SCORE = -1;


/*	static const int MATCH_SCORE = 1000;
	static const int MISMATCH_SCORE = -1000;
	static const int GAP_OPENING_SCORE = -1000;
	static const int GAP_EXTENSION_SCORE = -1;
*/
	string a;
	string b;
	string c;
	int **score;
	int **gapa;
	int **gapb;
	double identity;
	int p_start;
	int p_end;
	int a_start;
	int a_end;
	int anchor_len;
	int left_anchor;
	int right_anchor;

private:
	void clear(int, int);
	int max3(int, int, int);
	void print_matrix(string, const string &, const string &, int **);	

public:
	aligner(int anchor_length = 1,int reflen=10000);
	~aligner();
	void align(const string &, const string &);
	void extract_calls(const string&, int, int, int*, FILE*);
	void dump(FILE *fo);
	int get_start();
	int get_end();
	float get_identity();
	int get_left_anchor();
	int get_right_anchor();
	string get_a();
	string get_b();

};
#endif

