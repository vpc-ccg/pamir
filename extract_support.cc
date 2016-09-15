#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "common.h"
using namespace std;

void extract_support(string input_file, string out_prefix)
{
	char *line 			= new char[1000000];
	char *readName 		= new char[1000];
	char *flag 			= new char[10];
	char *contigName 	= new char[100];
	char *tmp 			= new char[10000];
	char *seq 			= new char [1000];
	char *rvSeq			= new char[1000];
	char *qual 			= new char[1000];
	char *rQual		= new char[1000];
	
	FILE *fin			= fopen(input_file.c_str(),"r");
	char *pContigName   = new char[100];
	strcpy(pContigName, "");
	fgets(line, 100000,fin);
	readName 			= strtok(line,"\t");	
	flag				= strtok(NULL,"\t");
	contigName 			= strtok(NULL,"\t");
	FILE *fout 			= fopen((out_prefix + "/" + string(contigName)+ "_reads.fq").c_str(),"w");
	FILE *fout2 		= fopen((out_prefix + "/" + string(contigName)+ "_rv_reads.fq").c_str(),"w");
	strcpy(pContigName,contigName);
	for(int i = 0; i<6; i++)
		tmp 			= strtok(NULL,"\t");
	seq 				= strtok(NULL,"\t");
	strcpy(rvSeq, reverse_complement(string(seq)).c_str());
	qual				= strtok(NULL,"\t");
	strcpy(rQual, reverse(string(qual)).c_str());
	fprintf(fout,"@%s\n%s\n+\n%s\n", readName, seq, qual);
	fprintf(fout2,"@%s\n%s\n+\n%s\n", readName, rvSeq, rQual);
	while(fgets(line, 100000,fin)!=NULL)
	{
		readName 		= strtok(line,"\t");	
		flag			= strtok(NULL,"\t");
		contigName 		= strtok(NULL,"\t");
		if(strcmp(pContigName,contigName)!=0)
		{
			fclose(fout);
			fclose(fout2);
			fout 		= fopen((out_prefix + "/" + string(contigName)+ "_reads.fq").c_str(),"w");
			fout2 		= fopen((out_prefix + "/" + string(contigName)+ "_rv_reads.fq").c_str(),"w");
			strcpy(pContigName,contigName);
		}
		for(int i = 0; i<6; i++)
			tmp 		= strtok(NULL,"\t");
		seq 			= strtok(NULL,"\t");
		strcpy(rvSeq, reverse_complement(string(seq)).c_str());
		qual			= strtok(NULL,"\t");
		strcpy(rQual, reverse(string(qual)).c_str());
		fprintf(fout,"@%s\n%s\n+\n%s\n", readName, seq, qual);
		fprintf(fout2,"@%s\n%s\n+\n%s\n", readName, rvSeq, rQual);
	}
	fclose(fin);
	fclose(fout);
	fclose(fout2);
}
int main(int argc,char *argv[])
{
	extract_support(string(argv[1]), string(argv[2]));
	return 0;
}
