#include <string>
#include "common.h"
#include "extractor.h"
#include "record.h"
#include "sam_parser.h"
#include "bam_parser.h"

using namespace std;


inline void output_record(FILE *fp, int ftype, const Record &rc)
{
	string record;
	uint32_t flag = rc.getMappingFlag();
	
	string mate = ((flag & 0x40)==0x40)?"/1":"/2";
	string seqname = rc.getReadName()+ mate;
	int reversed = ((flag & 0x10) == 0x10);
	string seq = (reversed) ? reverse_complement (rc.getSequence()) : rc.getSequence();
	string qual = (reversed) ? reverse (rc.getQuality()): rc.getQuality();

	if (ftype == 1)
		record = S(">%s\n%s\n", seqname.c_str(), seq.c_str());
	else if (ftype==2)
		record = S("@%s\n%s\n+\n%s\n", seqname.c_str(), seq.c_str(), qual.c_str());
	else if (ftype==3)
		record = S("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n",
				rc.getReadName(),
				rc.getMappingFlag(),
				rc.getChromosome(),
				rc.getLocation(),
				rc.getMappingQuality(),
				rc.getCigar(),
				rc.getPairChromosome(),
				rc.getPairLocation(),
				rc.getTemplateLength(),
				rc.getSequence(),
				rc.getQuality(),
				rc.getOptional()
				);

	fwrite(record.c_str(), 1, record.size(), fp);
}
/****************************************************************/
extractor::extractor(string filename, string output_prefix, int ftype, int oea, int orphan) 
{

	FILE *fi = fopen(filename.c_str(), "rb");
	char magic[2];
	fread(magic, 1, 2, fi);
	fclose(fi);

	string extensions[] = {"","fa","fq","sam"};

	Parser *parser;
	if (magic[0] == char(0x1f) && magic[1] == char(0x8b)) 
		parser = new BAMParser(filename);
	else
		parser = new SAMParser(filename);

	string comment = parser->readComment();


	FILE *foea_mapped, *foea_unmapped, *forphan, *fall_int;
	fall_int  = fopen ((output_prefix + "all_interleaved.fastq").c_str(), "w");
	// open file
	if (oea) 
	{
		string foea_mapped_name = output_prefix + "oea.mapped."+extensions[ftype];
		string foea_unmapped_name = output_prefix + "oea.unmapped."+extensions[ftype];
		foea_mapped  = fopen (foea_mapped_name.c_str(), "w");
		foea_unmapped  = fopen (foea_unmapped_name.c_str(), "w");
	}
		
	if (orphan)
	{
		string forphan_name = output_prefix + "orphan."+extensions[ftype];
		forphan = fopen (forphan_name.c_str(), "w");
	}

	
	uint32_t flag;
	while (parser->hasNext())
	{
		const Record &rc = parser->next();
		flag = rc.getMappingFlag();
		if ((flag & 0xD) == 0xD)
		{
			output_record (forphan, ftype, rc);
			output_record(fall_int, ftype, rc);
		}
		else if ((flag & 0x5) == 0x5) 
		{
			output_record (foea_unmapped, ftype, rc);
			output_record(fall_int, ftype, rc);
		}
		else if ((flag & 0x9) == 0x9)
		{
			output_record (foea_mapped, ftype, rc);
			output_record(fall_int, ftype, rc);
		}
		parser->readNext();
	}
	
	delete parser;
	
	// close file
	if (oea)
	{
		fclose(foea_mapped);
		fclose(foea_unmapped);
	}
	if (orphan)
	{
		fclose(forphan);
	}
	fclose(fall_int);
}
/***************************************************************/
int md_length( char *md)
{
	md+=5;
	int length = 0;
	int tmp = 0;
	while( *md )
	{
		if (isdigit(*md))
		{ 
			tmp = 10 * tmp + (*md - '0');
		}
		else
		{
			length += tmp;
			tmp = 0;
		}
		md++;
	}
	if (0 < tmp){length+=tmp;}
	return length;
}
/*************************************************************/
extractor::extractor(string filename, string output) 
{

	FILE *fi = fopen(filename.c_str(), "rb");
	char magic[2];
	fread(magic, 1, 2, fi);
	fclose(fi);

	Parser *parser;
	if (magic[0] == char(0x1f) && magic[1] == char(0x8b)) 
		parser = new BAMParser(filename);
	else
		parser = new SAMParser(filename);

	string comment = parser->readComment();


	FILE *fq = fopen (output.c_str(),"w");
	//FILE *fqsam = fopen ((output + ".sam").c_str(),"w");

	char *opt 			= new char[1000000];
	char *MD 			= new char[1000000];
	string fr_rname;
	string fr_seq;
	string fr_qual;
	uint32_t fr_flag;
	uint32_t flag;
	uint32_t tlen;
	int errNum;
	while (parser->hasNext())
	{
		const Record &rc = parser->next();
		flag = rc.getMappingFlag();
		tlen = strlen(rc.getSequence());
		errNum = int(tlen*0.94);
		if((flag & 0x800) != 0x800)
		{
			if((flag & 0x5) == 0x5 || (flag & 0x9)==0x9 || (flag & 0xD)==0xD)
			{
				output_record(fq, 2, rc);
			}
			else
			{
				strcpy(opt,rc.getOptional());
				strtok(opt,"\t");
				MD = strtok(NULL,"\t");
				if(md_length(MD) >= errNum)
				{
					if((flag & 0x40) == 0x40)
					{
						fr_rname = string(rc.getReadName());
						fr_seq = string(rc.getSequence());
						if((flag & 0x10) == 0x10)
						{
							fr_seq = reverse_complement(fr_seq);
						}
						fr_qual = string(rc.getQuality());
					}
				}
				else
				{
					if((flag & 0x40) == 0x40)
					{
						output_record(fq, 2, rc);
						//output_record(fqsam, 3, rc);
						int secondWritten = 0;
						while(secondWritten != 1)
						{
							parser->readNext();
							if(parser->hasNext())
							{
								const Record &rc2 = parser->next();
								flag = rc2.getMappingFlag();
								if((flag & 0x800) != 0x800)
								{
									output_record(fq, 2, rc2);
									//output_record(fqsam, 3, rc2);
									secondWritten = 1;
								}
							}
						}
					}
					else
					{
						fprintf(fq,"@%s/1\n%s\n+\n%s\n", fr_rname.c_str(), fr_seq.c_str(), fr_qual.c_str());
						//fprintf(fqsam,"@%s/1\t%s\t%s\n", fr_rname.c_str(), fr_seq.c_str(), fr_qual.c_str());
						output_record(fq, 2, rc);
						//output_record(fqsam, 3, rc);
					}
				}
			}
		}
		parser->readNext();
	}
	
	delete parser;
	fclose(fq);
	//fclose(fqsam);
}
/***************************************************************/
extractor::~extractor()
{
		
}
