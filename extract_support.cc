#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;

void extract_support(string input_file, string out_prefix)
{
	char *line 			= new char[1000000];
	char *readName 		= new char[1000];
	char *flag 			= new char[10];
	char *contigName 	= new char[100];
	char *tmp 			= new char[10000];
	char *seq 			= new char [1000];
	char *qual 			= new char[1000];
	
	FILE *fin			= fopen(input_file.c_str(),"r");
	char *pContigName   = new char[100];
	strcpy(pContigName, "");
	fgets(line, 100000,fin);
	readName 	= strtok(line,"\t");	
	flag		= strtok(NULL,"\t");
	contigName 	= strtok(NULL,"\t");
	FILE *fout = fopen((out_prefix + string(contigName)+ "_reads.fq").c_str(),"w");
	strcpy(pContigName,contigName);
	for(int i = 0; i<6; i++)
		tmp 		= strtok(NULL,"\t");
	seq 		= strtok(NULL,"\t");
	qual		= strtok(NULL,"\t");
	fprintf(fout,"@%s\n%s\n+\n%s\n", readName, seq, qual);
	while(fgets(line, 100000,fin)!=NULL)
	{
		readName 	= strtok(line,"\t");	
		flag		= strtok(NULL,"\t");
		contigName 	= strtok(NULL,"\t");
		if(strcmp(pContigName,contigName)!=0)
		{
			fclose(fout);
			fout = fopen((out_prefix + string(contigName)+ "_reads.fq").c_str(),"w");
			strcpy(pContigName,contigName);
		}
		for(int i = 0; i<6; i++)
			tmp 		= strtok(NULL,"\t");
		seq 		= strtok(NULL,"\t");
		qual		= strtok(NULL,"\t");
		fprintf(fout,"@%s\n%s\n+\n%s\n", readName, seq, qual);
	}
	fclose(fin);
	fclose(fout);
}
int main(int argc,char *argv[])
{
	extract_support(string(argv[1]), string(argv[2]));
	return 0;
}
