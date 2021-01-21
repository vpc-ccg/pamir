#include <string>
#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include "common.h"

using namespace std;

/**************************************************/
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
/**********************************************************/
string reverse_complement ( const string &str )
{
	const char *revComp = "TBGDEFCHIJKLMNOPQRSAUVWXYZ";

	string x;
	for (int i = str.size() - 1; i >= 0; i--)
		x += revComp[str[i] - 'A'];
	return x;
}
/*******************************************************************/
string reverse ( const string &str )
{
	string x;
	for (int i = str.size() - 1; i >= 0; i--)
		x += str[i];
	return x;
}
/*******************************************************************/
string S (const char* fmt, ...) {
	char *ptr = 0;
    va_list args;
    va_start(args, fmt);
    vasprintf(&ptr, fmt, args);
    va_end(args);
    string s = ptr;
    free(ptr);
    return s;
}
/********************************************************************/
int check_AT_GC(const string &contig, const double &MAX_AT_GC)
{
	double AT_count = 0, GC_count = 0, A_count = 0, T_count = 0, G_count = 0, C_count = 0;
	int clen        = contig.length();
	for(int i = 0; i < clen-1; i++)
	{
		if((contig[i] == 'A' && contig[i+1] == 'T') ||(contig[i] == 'T' && contig[i+1] == 'A')) AT_count++;
		if((contig[i] == 'C' && contig[i+1] == 'G') ||(contig[i] == 'G' && contig[i+1] == 'C'))	GC_count++;
		if((contig[i] == 'A' && contig[i+1] == 'A')) A_count++;
		if((contig[i] == 'T' && contig[i+1] == 'T')) T_count++;
		if((contig[i] == 'G' && contig[i+1] == 'G')) G_count++;
		if((contig[i] == 'C' && contig[i+1] == 'C')) C_count++;
	}
	if(AT_count/(double)clen >= MAX_AT_GC || GC_count/(double)clen >= MAX_AT_GC || A_count/(double)clen >=MAX_AT_GC || T_count/(double)clen >=MAX_AT_GC || G_count/(double)clen >=MAX_AT_GC || C_count/(double)clen >=MAX_AT_GC )
		return 0;
	else
		return 1;
}
/*******************************************************************/
void append_vcf(const string &chrName, const string &reference,
                const vector< tuple< string, int, int, string, int, float > > &reports, const int &clusterId,
                string &vcf_str, string &vcf_str_del ) {
    for (int r = 0; r < reports.size(); r++) {
        if(get<0>(reports[r]) == "INS") {
            vcf_str += 	chrName;	vcf_str += 	"\t";
            vcf_str +=	std::to_string(get<1>(reports[r]));	vcf_str += "\t.\t";
            vcf_str +=  reference.at(r);
            vcf_str +=  "\t";
            vcf_str += reference.at(r);
            vcf_str += get<3>(reports[r]);
            vcf_str +=  "\t";
            vcf_str +=  std::to_string( get<5>(reports[r]));
            vcf_str +=  "\tPASS\t";
            vcf_str += "Cluster=";	vcf_str +=	std::to_string(clusterId) ;
            vcf_str += ";Support=";	vcf_str	+=	std::to_string(get<4>(reports[r])) ;
            vcf_str += "\n";
        }
    }
}
/*******************************************************************/
void append_vcf_hybrid(const string &chrName, const string &reference,
                       const vector< tuple< string, int, int, string, int, float > > &reports, const int &clusterId, const string& cluster_type,
        int bimodals, int left, int right, int misc, int estimation,
        string &vcf_str, string &vcf_str_del ) {
    for (int r = 0; r < reports.size(); r++) {
        if(get<0>(reports[r]) == "INS") {
            vcf_str += 	chrName;	vcf_str += 	"\t";
            vcf_str +=	std::to_string(get<1>(reports[r]));	vcf_str += "\t.\t";
            vcf_str +=  reference.at(r);
            vcf_str +=  "\t";
            vcf_str += reference.at(r);
            vcf_str += get<3>(reports[r]);
            vcf_str +=  "\t";
            vcf_str +=  std::to_string( get<5>(reports[r]));
            vcf_str +=  "\tPASS\t";
            vcf_str += "Cluster=";	vcf_str +=	std::to_string(clusterId) ;
            vcf_str += ";Support=";	vcf_str	+=	std::to_string(get<4>(reports[r])) ;
            vcf_str += ";Type="; vcf_str += cluster_type;
            vcf_str += ";Bimodal="; vcf_str += std::to_string(bimodals);
            vcf_str += ";G1(L)="; vcf_str += std::to_string(left);
            vcf_str += ";G2(R)="; vcf_str += std::to_string(right);
            vcf_str += ";Misc="; vcf_str += std::to_string(misc);
            vcf_str += ";Estimate="; vcf_str += std::to_string(estimation);
            vcf_str += "\n";
        }
    }
}