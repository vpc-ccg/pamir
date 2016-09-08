#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <map>
using namespace std;

int main(int argc, char* argv[])
{
	if(argc!=3)
	{
	  printf("Usage: ./removeDuplicateInsertions inputFileName inputFileName.woDups)");
	  return 0;
	}
	FILE *fvcf;
	FILE *fout;
	fvcf=fopen(argv[1],"r");
	fout=fopen(argv[2],"w");
	char *str=new char[1000000]; 
	char *nextStr=new char[1000000]; 
	char *cpstr=new char[1000000]; 
	char *cpnextStr=new char[1000000]; 
	char *chrName1=new char[100];
	char *cpchrName1=new char[100];
	char *chrName2=new char[100];
	char *loc1=new char[100];
	char *loc2=new char[100];
	char *length1=new char[100];
	char *length2=new char[100];
	char *dot = new char[100];
	char *seq1=new char[100000];
	char *cpseq1=new char[100000];
	char *seq2=new char[100000];
	int iloc1;
	int iloc2;
	int ilength1;
	int ilength2;
	fgets(str, 1000000, fvcf);
	strcpy(cpstr,str);
	while(!feof(fvcf)&&fgets(nextStr, 1000000, fvcf)!=NULL)
	{	
		strcpy(cpnextStr,nextStr);
		chrName1=strtok(str,"\t");
		strcpy(cpchrName1,chrName1);
		loc1=strtok(NULL,"\t");
		iloc1=atoi(loc1);
		length1=strtok(NULL,"\t");
		ilength1=atoi(length1);
		dot = strtok(NULL,"\t");
		seq1=strtok(NULL,"\t");
		strcpy(cpseq1,seq1);
		chrName2=strtok(nextStr,"\t");
		loc2=strtok(NULL,"\t");
		iloc2=atoi(loc2);
		length2=strtok(NULL,"\t");
		ilength2=atoi(length1);
		dot =  strtok(NULL,"\t");
		seq2=strtok(NULL,"\t");
		if(strcmp(chrName1,chrName2)==0 && (iloc2==iloc1) && (ilength1==ilength2) && strcmp(cpseq1,seq2)==0)
		{
			strcpy(str,cpnextStr);
			strcpy(cpstr,str);
		}
		else
		{
			fprintf(fout,"%s\t%d\t%d\t%s\n",cpchrName1, iloc1, ilength1, cpseq1);
			strcpy(str,cpnextStr);
			strcpy(cpstr,str);
		}
	}
	char *tmp = new char[1000];
	tmp = strtok(cpstr, "\t");
	for(int i =0;i<3;i++)
	{
		fprintf(fout,"%s\t",tmp);
		tmp = strtok(NULL, "\t\n");
	}
	tmp = strtok(NULL, "\t\n");
	fprintf(fout,"%s\t",tmp);

	fclose(fvcf);
	fclose(fout);
	return 0;
}
