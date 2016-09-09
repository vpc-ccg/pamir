#include<iostream>
#include<string>
#include<algorithm>
#include<cmath>
#include<vector>
#include "aligner.h"

using namespace std;

aligner::aligner(int anchor_length,int reflen)
{
	MAX_SIDE=reflen;
	anchor_len = anchor_length;
	score = new int*[MAX_SIDE+1];
	gapa  = new int*[MAX_SIDE+1];
	gapb  = new int*[MAX_SIDE+1];

	for (int i=0; i<= MAX_SIDE; i++)
	{
		score[i] = new int[MAX_SIDE+1];
		gapa[i] = new int[MAX_SIDE+1];
		gapb[i] = new int[MAX_SIDE+1];
	}

	score[0][0] = gapa[0][0] = gapb[0][0] = 0;
	for (int i=1; i<= MAX_SIDE; i++)
	{
		score[i][0] = 0;
		gapa[0][i] = gapa[i][0] = GAP_OPENING_SCORE + (i-1)*GAP_EXTENSION_SCORE;
		gapb[0][i] = gapb[i][0] = GAP_OPENING_SCORE + (i-1)*GAP_EXTENSION_SCORE;

	}
	for (int i=1; i<=MAX_SIDE; i++)
		score[0][i] = GAP_OPENING_SCORE + (i-1)*GAP_EXTENSION_SCORE;

}

aligner::~aligner()
{
	for (int i=0; i<= MAX_SIDE; i++)
	{
		delete score[i];
		delete gapa[i];
		delete gapb[i];
	}	
	delete score;
}

void aligner::print_matrix(string name, const string &a, const string &b, int **matrix)
{
	fprintf(stdout, "=================================\n");
	fprintf(stdout, "Matrix: %s\n", name.c_str());
	fprintf(stdout, "=================================\n");
	fprintf(stdout,"       -");
	for (int i=0; i<a.length();i++)
		fprintf(stdout, "%4c", a[i]);
	fprintf(stdout, "\n");
	for (int j=0; j<=b.length(); j++)
	{
		if (j>0) fprintf(stdout,"%4c", b[j-1]);
		else fprintf(stdout, "   -");
		for (int i=0; i<=a.length();i++)
		{
			if (matrix[i][j] <= -100000)
				fprintf(stdout, "   N");
			else 
				fprintf(stdout,"%4d", matrix[i][j]);
		}
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "=================================\n");
}

void aligner::clear(int a, int b)
{
	identity = 0;

	for (int i=0; i<a; i++)
		for (int j=0; j<b; j++)
			score[i][j]=0;

	for (int i=0; i<b; i++)
	{	
		gapa[0][i]=-100000;
		for (int j=1;j<a;j++)
			gapa[j][i]=0;
	}
	for (int j=0; j<a; j++)
	{	
		gapb[j][0]=-100000;
		for (int i=1;i<b;i++)
			gapb[j][i]=0;
	}
}

int aligner::max3(int a, int b, int c)
{
	return max (max(a,b), c);
}

void aligner::align(const string &ref, const string &ass)
{
	int len = ass.length();

	int tmp;	
	for (int i=1; i<=ref.length();i++)
	{
		for (int j=1; j<= len;j++)
		{
			gapa[i][j] = max( gapa[i-1][j] + GAP_EXTENSION_SCORE, score[i-1][j]+GAP_OPENING_SCORE+GAP_EXTENSION_SCORE);
			gapb[i][j] = max( gapb[i][j-1] + GAP_EXTENSION_SCORE, score[i][j-1]+GAP_OPENING_SCORE+GAP_EXTENSION_SCORE);
			tmp = ((ref[i-1]==ass[j-1])?MATCH_SCORE:MISMATCH_SCORE);
			score[i][j] = max3(score[i-1][j-1]+tmp, gapa[i-1][j-1]+tmp, gapb[i-1][j-1]+tmp);
		}
	}
	
	//backtracking
	int max = score[1][len];
	int pi = 1;
	int pj = len;
	for (int i=1; i<=ref.length(); i++)
	{
		if (score[i][len] > max)
		{
			max = score[i][len];
			pi = i;
		}
	}

	p_end = pi-1;

	a_start = 0;
	a_end = pj-1;

	a = b = c  = "";

	int cur = 0;

	int match = 0, mismatch = 0, indel = 0;
	
	while ( pi >0 && pj > 0 )
	{
		if ( cur == 0 )
		{
			tmp = ((ref[pi-1] == ass[pj-1])?MATCH_SCORE:MISMATCH_SCORE);
			a = ref[pi-1] + a;
			b = ass[pj-1] + b;
			c = ((tmp>0)?'|':' ')+c;

			if (gapa[pi-1][pj-1] + tmp == score[pi][pj]) 
			{
				indel++;
				cur = 1;
			}
			else if (gapb[pi-1][pj-1] + tmp == score[pi][pj]) 
			{
				indel++;
				cur = 2;
			}
			pi--;
			pj--;

			if (tmp > 0){
				match++;
			}
			else{
				mismatch++;
			}
		}
		else if (cur == 1)
		{
			a = ref[pi-1] + a;
			b = '-' + b;
			c = ' '+c;
			if (score[pi-1][pj]+GAP_EXTENSION_SCORE+GAP_OPENING_SCORE == gapa[pi][pj])
				cur = 0;
			pi--;
		}
		else 
		{
			a = '-' + a;
			b = ass[pj-1] + b;
			c = ' '+c;
			if (score[pi][pj-1]+GAP_EXTENSION_SCORE+GAP_OPENING_SCORE == gapb[pi][pj])
				cur = 0;
			pj--;
		}
	}

/*	left_anchor = 0;
	for (int i=0; i<a.length(); i++)
		if (a[i] == b[i])
			left_anchor++;
		else
			break;
	right_anchor = 0;
	for (int i=a.length()-1; i>=0; i--)
		if (a[i] == b[i])
			right_anchor++;
		else
			break;
*/
	left_anchor = 0;
	right_anchor = 0;
	int mm =0;
	int totalErr=0;
	for (int i=0; i<a.length(); i++)
	{
		if (a[i] == b[i]){
			left_anchor++;
			mm = 0;
		}
		else
		{
			mm++;
			totalErr++;
			left_anchor++;
		}
		if (mm==2)
		{
			left_anchor-=2;
			break;
		}
		if( i > 10 && totalErr > (float(i)/5.0))
			break;
	}
	mm=0;
	totalErr=0;
	for (int i=a.length()-1; i>=0; i--)
	{
		if (a[i] == b[i]){
			right_anchor++;
			mm = 0;
		}
		else
		{
			mm++;
			totalErr++;
			right_anchor++;
		}
		if(mm==2)
		{
			right_anchor-=2;
			break;
		}
		if(i < a.length()-10 && totalErr > (float(a.length()-float(i))/5.0))
			break;
	}

	p_start = pi;
	len = ass.length();
	int l = 1+abs((((p_end-p_start)+1)-len));
	identity = 1 - (mismatch+indel+log(l))/len;
}
/**********************************************************************/
int aligner::extract_calls( const int &cluster_id, vector<tuple<string, int, int, string, int, float > > &reports, const int &contig_support, const int &ref_start)
{
	dump(stdout);
	int mapped					= 0;
	int insertion_start_loc 	= -1;
	int insertion_end_loc		= -1;
	double fwdIden 				= get_identity();
	string insertion_content;
	int fwdS, fwdE, fwdAS, fwdAE;
	if( fwdIden > 1.0 && ( left_anchor > 5 && right_anchor > 5 ) && ( left_anchor > 16 || right_anchor > 16 ) )
	{
		fwdS = ref_start + p_start;
		fwdAS 					= 0;
		fwdAE 					= a.length()-1;
		int del_len				= 0;
		int ins_len				= 0;
		int deletion_start_loc 	= 0;
		int deletion_end_loc 	= 0;
		int isdeletion			= 0;
		string deletion_content;
		cout<<"fwdIden: "<<fwdIden<<"\nfwdAE: "<<fwdAE<<"\tfwdAS: "<<"\tfwdS: "<<fwdS<<"\t"<<p_start<<endl;			
		cout<<"REF PART OF MAPPING: "<<a.c_str()<<endl;
		if( b.length() > 0 )
		{
			int p = 0;
			while( p < b.length() )
			{
				del_len=0;
				for( ; p < b.length(); p++ )
				{
					if( b[p] == '-' )
					{
						deletion_start_loc = p;
						break;
						del_len = 1;
					}
				}
				for( ; p < b.length(); p++ )
				{
					if( b[p] != '-' )
					{
						deletion_end_loc = p;
						break;
					}
					else
						del_len++;
				}
				if( del_len >= 5 && deletion_start_loc != -1 && deletion_start_loc + 1 != insertion_end_loc )
				{
					isdeletion         = 1;
					deletion_content   = a.substr(deletion_start_loc,deletion_end_loc-deletion_start_loc);
					deletion_start_loc = deletion_start_loc+fwdS;
					if(deletion_content.length()>0)
					{
						reports.push_back(tuple<string, int, int, string, int, float>("DEL", deletion_start_loc, deletion_content.length(), deletion_content, contig_support, identity ) );
					}
				}
			}
		}
		if( isdeletion == 0 and a.length() > 0 )
		{
			int p 						= 0;
			int decreaseFromInsStartLoc = 0;
			while( p < a.length() )
			{
				ins_len = 0;
				for( ; p < a.length(); p++ )
				{
					if( a[p] == '-' )
					{
						insertion_start_loc = p;
						break;
						ins_len = 1;
					}
				}
				for( ; p < a.length(); p++ )
				{
					if( a[p] != '-' )
					{
						insertion_end_loc = p;
						break;
					}
					else
						ins_len++;
				}
				if( ins_len >= 5 && insertion_start_loc != -1 && insertion_start_loc + 1 != insertion_end_loc )
				{
					insertion_content   = b.substr(insertion_start_loc,insertion_end_loc-insertion_start_loc);
					insertion_start_loc = insertion_start_loc+fwdS;
					if( insertion_content.length() > 0 )
					{
						cout<<insertion_start_loc<<"\t"<<insertion_content.length()<<"\t"<<insertion_content.c_str()<<"\t"<<contig_support<<endl;
						cout<<insertion_content.length()<<endl;
						reports.push_back(tuple<string, int, int, string, int, float>("INS", insertion_start_loc, insertion_content.length(), insertion_content, contig_support, identity ) );
						mapped = 1;
					}
				}
				p++;
			}
		}
	}
	return mapped;
}
/********************************************************************************/
void aligner::dump(FILE *fo)
{
	string cnt = "";
	for (int i=1; i<=a.length(); i++)
	{
		if (i%10==0)
			cnt+='0';
		else if(i%5==0)
			cnt+='5';
		else
			cnt+=' ';
	}
	fprintf(fo, "               %s\nG(%10d): %s\n               %s\nA(%10d): %s\n",cnt.c_str(), p_start,a.c_str(), c.c_str(), 1,b.c_str());
	//fprintf(fo, "   %s\nG: %s\n   %s\nA: %s\n",cnt.c_str(), a.c_str(), c.c_str(), b.c_str()); 
}
/********************************************************************************/
int aligner::get_start()
{
	return p_start;
}
/********************************************************************************/
int aligner::get_end()
{
	return p_end;
}
/********************************************************************************/
float aligner::get_identity()
{
	return 100*identity;
}
/********************************************************************************/
string aligner::get_a()
{
	return a;
}
/********************************************************************************/
string aligner::get_b()
{
	return b;
}
/********************************************************************************/
int aligner::get_left_anchor()
{
	return left_anchor;
}
/********************************************************************************/
int aligner::get_right_anchor()
{
	return right_anchor;
}
