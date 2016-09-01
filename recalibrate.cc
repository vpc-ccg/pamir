#include <iostream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;

int b_search(vector<int> &starts, int t, int b, int e){
	if (b==e) {
		if (starts[b]<=t)
			return b;
		else 
			return b-1;
	}
	else {
		int m = (b+e)/2;
		if (t>=starts[m])
			return b_search(starts, t, m+1, e);
		else {
			return b_search(starts, t, b, m);
		}
	}

}

int main(int argc, char** argv) {

	if(argc!=4)
	{
		fprintf(stderr,"The number of inputs are not correct! Check your command!\n");
		exit(1);
	}
	vector<int> starts;
	vector<string> names;
	ifstream fin;
	int start;
	string name;

	string fname = argv[1];
	fin.open(fname.c_str());

	while (fin >> name >> start )
	{
		starts.push_back (start);
		names.push_back(name);
	}
	FILE *fin2 = fopen(argv[2],"r");
	FILE *fout = fopen(argv[3],"w");
	char *funmapped=new char[1000];
	sprintf(funmapped,"%s.unmapped",argv[3]);
	FILE *fout2 = fopen(funmapped,"w");
	char *line = new char[1000000];
	char *cpline = new char[1000000];
	char *readName = new char[100];
	char *flag     = new char[100];
	char *chrName  = new char[100];
	char *chrLoc   = new char[100];
	char *tmp      = new char[100];
	int ichrLoc;
	while (fgets(line,1000000,fin2)!=NULL)
	{
		readName=strtok(line,"\t");
		flag=strtok(NULL,"\t");
		chrName=strtok(NULL, "\t");
		chrLoc=strtok(NULL, "\t");
		ichrLoc=atoi(chrLoc);
		if(strcmp(chrName,"*")!=0)
		{
			int x =  b_search(starts, ichrLoc, 0, starts.size()-1);
			fprintf(fout,"%s\t%s\t%s\t%d",readName, flag, names[x].c_str(),ichrLoc-starts[x]+1);
			tmp=strtok(NULL,"\t\n");
			while(tmp!=NULL and strcmp(tmp,"")!=0)
			{
				fprintf(fout,"\t%s",tmp);
				tmp=strtok(NULL,"\t\n");
			}
			fprintf(fout,"\n");

		}
		else
		{
			fprintf(fout2,"%s\t%s\t*\t0",readName, flag);
			tmp=strtok(NULL,"\t\n");
			while(tmp!=NULL and strcmp(tmp,"")!=0)
			{
				fprintf(fout2,"\t%s",tmp);
				tmp=strtok(NULL,"\t\n");
			}
			fprintf(fout2,"\n");
		}
		delete line;
		line = new char[1000000];
	}
	fin.close();
	fclose(fin2);
	fclose(fout);
	fclose(fout2);
}
