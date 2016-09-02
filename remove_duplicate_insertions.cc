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
	char *mode1 = new char[100];
	char *mode2 = new char[100];
	char *chrName1=new char[100];
	char *cpchrName1=new char[100];
	char *chrName2=new char[100];
	char *loc1=new char[100];
	char *loc2=new char[100];
	char *length1=new char[100];
	char *length2=new char[100];
	char *seq1=new char[100000];
	char *cpseq1=new char[100000];
	char *seq2=new char[100000];
	char *support1=new char[100];
	char *support2=new char[100];
	int iloc1;
	int iloc2;
	int ilength1;
	int ilength2;
	int isupport1;
	int isupport2;
	fgets(str, 1000000, fvcf);
	strcpy(cpstr,str);
	int totalSupport=0;
	while(!feof(fvcf)&&fgets(nextStr, 1000000, fvcf)!=NULL)
	{	
		strcpy(cpnextStr,nextStr);
		mode1 = strtok(str,"\t");
		chrName1=strtok(NULL,"\t");
		strcpy(cpchrName1,chrName1);
		loc1=strtok(NULL,"\t");
		iloc1=atoi(loc1);
		length1=strtok(NULL,"\t");
		ilength1=atoi(length1);
		seq1=strtok(NULL,"\t");
		strcpy(cpseq1,seq1);
		support1=strtok(NULL,"\t\n");
		isupport1=atoi(support1);
		totalSupport+=isupport1;
		mode2 = strtok(nextStr,"\t");
		chrName2=strtok(NULL,"\t");
		loc2=strtok(NULL,"\t");
		iloc2=atoi(loc2);
		length2=strtok(NULL,"\t");
		ilength2=atoi(length1);
		seq2=strtok(NULL,"\t");
		support2=strtok(NULL,"\t\n");
		isupport2=atoi(support1);
		if(strcmp(chrName1,chrName2)==0 && (iloc2==iloc1) && (ilength1==ilength2) && strcmp(cpseq1,seq2)==0)
		{
			strcpy(str,cpnextStr);
			strcpy(cpstr,str);
		}
		else
		{
			fprintf(fout,"%s\t%d\t%d\t%s\t%d\n",cpchrName1, iloc1, ilength1, cpseq1, totalSupport);
			totalSupport=0;
			strcpy(str,cpnextStr);
			strcpy(cpstr,str);
		}
	}
	char *tmp = new char[1000];
	tmp = strtok(cpstr, "\t");
	tmp = strtok(NULL, "\t");
	while(tmp!=NULL)
	{
		fprintf(fout,"%s\t",tmp);
		tmp = strtok(NULL, "\t\n");
	}
	fclose(fvcf);
	fclose(fout);
	return 0;
}
