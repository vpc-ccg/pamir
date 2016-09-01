#include<iostream>
#include<string>
#include<algorithm>
#include<cmath>
#include"aligner.h"

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
	int pi=1;
	int pj=len;
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

	int match=0, mismatch=0, indel=0;
	
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
	for (int i=0; i<a.length(); i++)
	{
		if (a[i] == b[i]){
			left_anchor++;
			mm = 0;
		}
		else
		{
			mm++;
			left_anchor++;
		}
		if (mm==2)
		{
			left_anchor-=2;
			break;
		}
	}
	mm=0;
	for (int i=a.length()-1; i>=0; i--)
	{
		if (a[i] == b[i]){
			right_anchor++;
			mm = 0;
		}
		else
		{
			mm++;
			right_anchor++;
		}
		if(mm==2)
		{
			right_anchor-=2;
			break;
		}
	}

	p_start = pi;
	len = ass.length();
	int l = 1+abs((((p_end-p_start)+1)-len));
	identity = 1 - (mismatch+indel+log(l))/len;
}

void aligner::extract_calls(const string& ref, int start, int clusterId, int *coverage, FILE *fo)
{
	int g_pos = start;
	int r_pos = 1;
	int ev_gen_start, ev_gen_end, ev_ass_start, ev_ass_end;
	int i=0;
	string x, y, type;
	while (identity > .900 && i<c.length())
	{
		if (c[i]=='|')
		{
			r_pos++;
			g_pos++;
			i++;
		}
		else 
		{
			if (a[i] != '-' && b[i] != '-')
			{
				ev_gen_start = ev_gen_end = g_pos;
				ev_ass_start = ev_ass_end = r_pos;
				r_pos++;
				g_pos++;
				x=a[i];
				y=b[i];
				i++;
			}
			else if (a[i] == '-')
			{

				ev_gen_start = g_pos-1;
				ev_ass_start = r_pos-1;
				x=a[i-1];
				y=b[i-1];

				while (a[i]=='-')
				{
					y+=b[i];
					r_pos++;
					i++;
				}

				ev_gen_end = ev_gen_start;
				ev_ass_end=r_pos-1;
			}
			else if (b[i]=='-')
			{
				ev_gen_start = g_pos-1;
				ev_ass_start = r_pos-1;
				x = a[i-1];
				y = b[i-1];

				while (b[i]=='-')
				{
					x+=a[i];
					g_pos++;
					i++;
				}
				ev_gen_end = g_pos-1;
				ev_ass_end = ev_ass_start;
			}
			int min = 10000, max = -1;
			for (int kk=ev_ass_start-1; kk<ev_ass_end; kk++)
			{
				if (coverage[kk]>max)
					max = coverage[kk];
				if (coverage[kk]<min)
					min = coverage[kk];
			}

			if(x.length() > y.length()){
				type = "DEL";
			}
			else if(x.length() < y.length()){
				type = "INS";
			}
			else{
				type = "SNV";
			}

			fprintf(fo, "%s\t%d\t.\t%s\t%s\t%f\tPASS\tSVTYPE=%s;END=%d;Cluster=%d;MinSupport=%d;MaxSupport=%d;Identity=%f\n", 
					ref.c_str(), ev_gen_start, x.c_str(), y.c_str(), (-10*log(1-identity)), type.c_str(), ev_gen_end, clusterId, min, max, identity*100); 
		}
	}
}

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
	fprintf(fo, "   %s\nG: %s\n   %s\nA: %s\n",cnt.c_str(), a.c_str(), c.c_str(), b.c_str()); 
}

int aligner::get_start()
{
	return p_start;
}

int aligner::get_end()
{
	return p_end;
}

float aligner::get_identity()
{
	return 100*identity;
}

string aligner::get_a()
{
	return a;
}

string aligner::get_b()
{
	return b;
}

int aligner::get_left_anchor()
{
	return left_anchor;
}

int aligner::get_right_anchor()
{
	return right_anchor;
}
