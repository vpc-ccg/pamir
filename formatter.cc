#include<string>
#include "formatter.h"
formatter::formatter(string filename, int fqi, string output, int cnt, int fq) 
{
	// only on no interleaved files

	FILE *fp = fopen (filename.c_str(), "r");
	FILE *op = fopen (output.c_str(), "w");


	if (fqi == 0 )
		fq = 0;

	char st[500000];
	int size;
	char name[1000];
	while (NULL!=fgets(st, 500000, fp)) 
	{
		size = strlen(st);
		if (size > 2 && st[size - 3] == '/')
			fprintf(op, "%c%d/%c\n", (fq?'@':'>'), cnt, st[size-2]);
		else
			fprintf(op, "%c%d\n", (fq?'@':'>'), cnt);
		cnt++;


		fgets(st,500000,fp); 
		fprintf(op,"%s", st);
		if (fqi) {
			fgets(st,500000,fp); 
			if (fq)	fprintf(op,"%s", st);
			fgets(st,500000,fp); 
			if (fq) fprintf(op,"%s", st);
		}
	}
	fclose(fp);
	fclose(op);
}
