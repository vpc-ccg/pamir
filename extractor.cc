#include <string>
#include <DeeZ.h>
#include "extractor.h"

using namespace std;


inline string reverse_complete (const string &str)
{
	const char *revComp = "TBGDEFCHIJKLMNOPQRSAUVWXYZ";

	string x;
	for (int i = str.size() - 1; i >= 0; i--)
		x += revComp[str[i] - 'A'];
	return x;
}

inline string reverse (const string &str)
{
	string x;
	for (int i = str.size() - 1; i >= 0; i--)
		x += str[i];
	return x;
}


inline void output_record(FILE *fp, int ftype, const Record &rc)
{
	string record;
	uint32_t flag = rc.getMappingFlag();

	
	//string mate = ((flag & 0x40)==0x40)?"/1":"/2";
	string mate = ((flag & 0x40)==0x40)?"+":"-";
	string seqname = rc.getReadName()+ mate;
	int reversed = ((flag & 0x10) == 0x10);
	string seq = (reversed) ? reverse_complete (rc.getSequence()) : rc.getSequence();
	string qual = (reversed) ? reverse (rc.getQuality()): rc.getQuality();

	if (ftype == 1)
		record = S(">%s\n%s\n", seqname.c_str(), seq.c_str());
	else if (ftype==2)
		record = S("@%s\n%s\n+\n%s\n", seqname.c_str(), seq.c_str(), qual.c_str());
	else if (ftype==3)
		record = S("%s %d %s %d %d %s %s %d %d %s %s %s\n",
				rc.getReadName(),
				rc.getMappingFlag(),
				rc.getChromosome(),
				rc.getLocation(),
				rc.getMappingQuality(),
				rc.getCigar(),
				rc.getPairChromosome(),
				rc.getPairLocation(),
				rc.getTemplateLenght(),
				rc.getSequence(),
				rc.getQuality(),
				rc.getOptional()
				);

	fwrite(record.c_str(), 1, record.size(), fp);
}


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


	FILE *foea_mapped, *foea_unmapped, *forphan;
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
//		printf("%s\t%d %d\n", rc.getReadName(), flag, (flag &0x5));
		if ((flag & 0xD) == 0xD)
//			fprintf(forphan, "1\n");
			output_record (forphan, ftype, rc);
		else if ((flag & 0x5) == 0x5) 
//			fprintf(foea_unmapped, "1\n");
			output_record (foea_unmapped, ftype, rc);
		else if ((flag & 0x9) == 0x9)
//			fprintf(foea_mapped, "1\n");
			output_record (foea_mapped, ftype, rc);
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
}

extractor::~extractor()
{
		
}
