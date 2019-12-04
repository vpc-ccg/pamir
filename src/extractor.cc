#include <string>
#include <map>
#include <zlib.h>
#include <algorithm>
#include <map>
#include "common.h"
#include "extractor.h"
#include "record.h"
#include "sam_parser.h"
#include "stats.h"

using namespace std;
/****************************************************************/
inline void output_record(FILE *fp, int ftype, const Record &rc)
{
	string record;
	uint32_t flag = rc.getMappingFlag();
	
	string mate = ((flag & 0x40)==0x40)?"/1":"/2";
	int reversed = ((flag & 0x10) == 0x10);
	string seq = (reversed) ? reverse_complement (rc.getSequence()) : rc.getSequence();
	string qual = (reversed) ? reverse (rc.getQuality()): rc.getQuality();

	if (ftype == 1)
		record = S(">%s%s\n%s\n", rc.getReadName(), mate.c_str(), seq.c_str());
	else if (ftype==2)
		record = S("@%s%s\n%s\n+\n%s\n", rc.getReadName(), mate.c_str(), seq.c_str(), qual.c_str());
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
	else if (ftype==4)
		record = S("@%s\n%s\n+\n%s\n", rc.getReadName(),  seq.c_str(), qual.c_str());

	fwrite(record.c_str(), 1, record.size(), fp);
}
/****************************************************************/
inline void gz_record( gzFile fp, int ftype, const Record &rc)
{
	string record;
	uint32_t flag = rc.getMappingFlag();
	
	string mate = ((flag & 0x40)==0x40)?"/1":"/2";
	int reversed = ((flag & 0x10) == 0x10);
	string seq = (reversed) ? reverse_complement (rc.getSequence()) : rc.getSequence();
	string qual = (reversed) ? reverse (rc.getQuality()): rc.getQuality();

	if (ftype == 1)
		record = S(">%s%s\n%s\n", rc.getReadName(), mate.c_str(), seq.c_str());
	else if (ftype==2)
		record = S("@%s%s\n%s\n+\n%s\n", rc.getReadName(), mate.c_str(), seq.c_str(), qual.c_str());
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
	else if (ftype==4)
		record = S("@%s\n%s\n+\n%s\n", rc.getReadName(),  seq.c_str(), qual.c_str());

	//fwrite(record.c_str(), 1, record.size(), fp);
	gzwrite( fp, record.c_str(), record.size() );
}
/****************************************************************/
int parse_sc( const char *cigar, int &match_l, int &read_l )
{	
	int tmp = 0;
	match_l = 0; read_l = 0;

	while( *cigar )
	{
		if (isdigit(*cigar))
		{ tmp = 10 * tmp + (*cigar - '0');}
		else
		{
			if ( 'M' == *cigar )
			{	match_l += tmp;	}
			if ( 'D' != *cigar )
			{	read_l += tmp;	}
			else
			{ match_l -= tmp; }
			tmp = 0;
		}
		cigar++;
	}
	if (0>match_l){match_l=0;}
	return 0;
}

void output_mates(const Record &mate1, const Record &mate2, FILE *fout, int ftype) {
    u_int32_t FIRST_MATE = 0x40;
    if ((mate1.getMappingFlag() & FIRST_MATE) == FIRST_MATE) {
        output_record(fout, ftype, mate1);
        output_record(fout, ftype, mate2);
    } else {
        output_record(fout, ftype, mate2);
        output_record(fout, ftype, mate1);
    }
}

int is_good_mapping (const Record &mapping, float clip_ratio) {
    if ((mapping.getMappingFlag() & 0x4)== 0x4)
        return 0;

    int matched = 0;
    int length = 0;
    parse_sc( mapping.getCigar(), matched, length);

    if ( ( (float)matched/mapping.getSequenceSize() ) <= clip_ratio )
        return 0;
    else
        return 1;
}

extractor::extractor(string filename, string output_prefix, int ftype, int oea, int orphan, double clip_ratio = 0.99) {
    SAMParser *parser = new SAMParser(filename);
    string comment = parser->readComment();


    map<int, int> dist;
    map<int, int> read_lengths;
    map <string, Record> map_read;
    map <string, Record> map_orphan;
    map <string, Record> map_oea;

    int supp_cnt = 0;
    int orig_concordant = 0;
    int orig_discordant = 0;
    int orig_chimeric = 0;
    int orig_oea = 0;
    int orig_orphan = 0;

    int proc_concordant = 0;
    int proc_discordant = 0;
    int proc_chimeric = 0;
    int proc_oea = 0;
    int proc_orphan = 0;

    string extensions[] = {"", "fa", "fq", "sam"};

    FILE *foea_mapped, *foea_unmapped, *forphan, *forphan_sup, *fall_int, *f_min_length;
    fall_int = fopen((output_prefix + ".all_interleaved.fastq").c_str(), "w");
    if (oea) {
        string foea_mapped_name = output_prefix + ".oea.mapped." + extensions[ftype];
        string foea_unmapped_name = output_prefix + ".oea.unmapped." + extensions[ftype];
        foea_mapped = fopen(foea_mapped_name.c_str(), "w");
        foea_unmapped = fopen(foea_unmapped_name.c_str(), "w");
    }

    if (orphan) {
        string forphan_name = output_prefix + ".orphan.canonical." + extensions[ftype];
        string forphan_name_sup = output_prefix + ".orphan.almost." + extensions[ftype];
        forphan = fopen(forphan_name.c_str(), "w");
        forphan_sup = fopen(forphan_name_sup.c_str(), "w");
    }

    int count = 0;
    while (parser->hasNext()) {
        const Record &rc = parser->next();
        Record cur_rc(rc);


        string name = string(rc.getReadName());
        int flag = rc.getMappingFlag();
        int tlen = rc.getTemplateLength();

        if (flag < 256) // To-Do: include supplementary split-mapping as potential mapping locations
        {
            int length = rc.getSequenceSize();
            read_lengths[ length]++;

            if ((flag & 0xc) == 0xc) {
                // I am orphan
                orig_orphan++;
                proc_orphan++;

                auto it = map_orphan.find(name);
                if (it != map_orphan.end()) {
                    output_mates(rc, it->second, forphan, ftype);
                    output_mates(rc, it->second, fall_int, ftype);
                    map_orphan.erase(it);
                } else {
                    map_orphan[name] = rc;
                }
            } else if ( (flag & 0xc) == 0x4 || (flag & 0xc) == 0x8) {
                // I am OEA
                orig_oea++;

                auto it = map_oea.find(name);
                if (it != map_oea.end()) {
                    int gm1 = is_good_mapping(rc, clip_ratio);
                    int gm2 = is_good_mapping(it->second, clip_ratio);
                    if ( gm1+gm2 == 0) {
                        proc_orphan+=2;
                        output_mates(rc, it->second, forphan_sup, ftype);
                        output_mates(rc, it->second, fall_int, ftype);
                    } else {
                        proc_oea+=2;
                        if (gm1) {
                            output_record(foea_mapped, ftype, rc);
                            output_record(foea_unmapped, ftype, it->second);
                            output_record(fall_int, ftype, rc);
                            output_record(fall_int, ftype, it->second);
                        } else {
                            output_record(foea_unmapped, ftype, rc);
                            output_record(foea_mapped, ftype, it->second);
                            output_record(fall_int, ftype, rc);
                            output_record(fall_int, ftype, it->second);
                        }
                    }
                    map_oea.erase(it);
                } else {
                    map_oea[name] = rc;
                }
            } else {
                // I ROCK
                if ((flag & 0x2) == 2)
                    orig_concordant++;
                else {
                    if (strcmp(rc.getChromosome(), rc.getPairChromosome()) != 0 && strcmp(rc.getPairChromosome(), "=") !=0 )
                        orig_chimeric++;
                    else
                        orig_discordant++;
                }

                auto it = map_read.find(name);

                if (it != map_read.end()) {

                    if (strcmp(rc.getChromosome(), (it->second).getChromosome() ) != 0)
                    {
                        proc_orphan +=2;
                        output_mates(rc, it->second, forphan, ftype);
                        output_mates(rc, it->second, fall_int, ftype);
                    } else {
                        int gm1 = is_good_mapping(rc, clip_ratio);
                        int gm2 = is_good_mapping(it->second, clip_ratio);
                        if (gm1 + gm2 == 2) {
                            if ((flag & 0x2) == 0x2) {
                                dist[abs(tlen)]++;
                                proc_concordant+=2;
                            } else { // discordant
                                proc_discordant+=2;
// TODO: REMOVE this after tests.
                                output_record(foea_mapped, ftype, rc);
                                output_record(foea_unmapped, ftype, it->second);
                                output_record(fall_int, ftype, rc);
                                output_record(fall_int, ftype, it->second);
                            }
                        } else if (gm1 + gm2 == 0) {
                            proc_orphan+=2;
                            output_mates(rc, it->second, forphan_sup, ftype);
                            output_mates(rc, it->second, fall_int, ftype);
                        } else {
                            proc_oea+=2;
                            if (gm1) {
                                output_record(foea_mapped, ftype, rc);
                                output_record(foea_unmapped, ftype, it->second);
                                output_record(fall_int, ftype, rc);
                                output_record(fall_int, ftype, it->second);
                            } else {
                                output_record(foea_unmapped, ftype, rc);
                                output_record(foea_mapped, ftype, it->second);
                                output_record(fall_int, ftype, rc);
                                output_record(fall_int, ftype, it->second);
                            }
                        }
                    }
                    map_read.erase(it);
                } else {
                    map_read[name] = rc;
                }
            }
            count++;
            if (0 == count % 1000000) {
                fprintf(stderr, "\rProcessed %dM records.", count / 1000000);
            }
        } else {
            supp_cnt++;
        }
        parser->readNext();
    }
    fprintf(stderr, "\n");

    delete parser;

    // sanity checks
    if (0 < map_orphan.size()) {
        Logger::instance().error("Orphan Queue is not cleared:\n");
        map<string, Record>::iterator it;
        for (it = map_orphan.begin(); it != map_orphan.end(); it++) {
            Logger::instance().error("%s\n", it->second.getReadName());
        }
    }

    if (0 < map_read.size()) {
        Logger::instance().error("Big Queue  is not cleared:\n");
        map<string, Record>::iterator it;
        for (it = map_read.begin(); it != map_read.end(); it++) {
            Logger::instance().error("%s\n", it->second.getReadName());
        }
    }

    FILE *f_stat;
    f_stat = fopen((output_prefix + ".stat").c_str(), "w");


    fprintf(f_stat,"Total Number of Reads: %d\n",count / 2);
    fprintf(f_stat,"# of Primary Mappings: %d\n",count / 2);
    fprintf(f_stat,"  # of Supp. Mappings: %d\n\n",supp_cnt / 2);
    fprintf(f_stat,"Original:\n");
    fprintf(f_stat,"Concordant: %d\n",orig_concordant/2);
    fprintf(f_stat,"Discordant: %d\n",orig_discordant/2);
    fprintf(f_stat,"  Chimeric: %d\n",orig_chimeric/2);
    fprintf(f_stat,"       OEA: %d\n",orig_oea/2);
    fprintf(f_stat,"    Orphan: %d\n\n",orig_orphan/2);
    fprintf(f_stat,"Processed:\n");
    fprintf(f_stat,"Concordant: %d\n",proc_concordant/2);
    fprintf(f_stat,"Discordant: %d\n",proc_discordant/2);
    fprintf(f_stat,"  Chimeric: %d\n",proc_chimeric/2);
    fprintf(f_stat,"       OEA: %d\n",proc_oea/2);
    fprintf(f_stat,"    Orphan: %d\n\n",proc_orphan/2);
    if (read_lengths.size() > 0)
        fprintf(f_stat,"Read Length Range: [%d, %d]\n\n", read_lengths.begin()->first, read_lengths.rbegin()->first);
    if (dist.size() > 0) {
        fprintf(f_stat,"TLEN:\n");
        fprintf(f_stat,"Range: [%d, %d]\n", dist.begin()->first, dist.rbegin()->first);
        auto dist_stats = get_distribution(dist.begin(), dist.end(),
                                           [](const decltype(dist.begin()) &val) -> int { return val->first; });
        fprintf(f_stat," Mean: %.2lf\n", dist_stats.first);
        fprintf(f_stat,"  Std: %.2lf\n\n", dist_stats.second);
    }

    fprintf(f_stat,"Template Length Counts:\n");

    if (dist.size() > 0) {
        for (auto it=dist.begin(); it != dist.end(); it++)
            fprintf(f_stat,"%d\t%d\n", it->first, it->second);
    }


    fclose(f_stat);

    // close file
    if (oea) {
        fclose(foea_mapped);
        fclose(foea_unmapped);
    }
    if (orphan) {

        fclose(forphan);
        fclose(forphan_sup);
    }
    fclose(fall_int);
}
/***************************************************************/
extractor::~extractor()
{
}
