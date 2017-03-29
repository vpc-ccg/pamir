#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
using namespace std;

int main(int argc, char* argv[])
{
	if(argc!=3)
	{
	  printf("Usage: ./clean megablastFile cleanFile\n");
	  exit(0);
	}
        FILE *fin;
        FILE *fout;
	FILE *fout2;
	FILE *fout3;
	char *line =new char[1000000];
	char *line1 =new char[1000000];
	char *line2 =new char[1000000];
	char *line3 =new char[1000000];
	char *line4 =new char[1000000];
	/*********DELETE NEW LINES***********/
	fout2=fopen(argv[2],"w");
	fin=fopen(argv[1],"r");
	while(!feof(fin))
	{
		fgets (line, 1000000, fin);
		if((line[0]==' ' && line[1]==' ' && line[2]=='D' && line[3]=='a' && line[4]=='t' && line[5]=='a' && line[6]=='b' && line[7]=='a'))
		{
			break;
		}
		while(!(line[0]=='Q' && line[1]=='u' && line[2]=='e' && line[3]=='r' && line[4]=='y' && line[5]=='='))
		{	      
		       fgets (line, 1000000, fin);
		}
		fprintf(fout2,"%s",line);
		fgets (line, 1000000, fin);
		fgets (line, 1000000, fin);
		fprintf(fout2,"%s",line);
		fgets (line, 1000000, fin);
		fgets (line, 1000000, fin);
		if(line[0]=='S' && line[1]=='e' && line[2]=='q')
		{
			fgets(line, 1000000, fin);
			int a=0;
			while(line[0]!='>')
			{
				a++;
				fgets(line, 1000000, fin);
				if(a<=5 && strcmp(line,"\n")!=0 && line[0]!='>')
					fprintf(fout2,"%s",line);
			}
		}
		while(!(line[0]=='E' && line[1]=='f' && line[2]=='f' && line[3]=='e' && line[4]=='c' && line[5]=='t'))
		{
			fgets (line, 1000000, fin);
		}
		fgets (line, 1000000, fin);
		fgets (line, 1000000, fin);		
	}
        fclose(fin);
	fclose(fout2);
	return 0;
}
