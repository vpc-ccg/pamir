#ifndef SAM_PROCCESING_H
#define SAM_PROCCESING_H

#include "cigar.h"
#include <array>
#include <string>
#include <map>
#include <cstring>
enum class read_type{
    bpsup,
    bpsupoea,
    bpsuporphan,
    refsup,
    falsesup,
    nosup,
    unset,
};


std::string &get_read_type_str(read_type type);

enum class cigar_character_type{
    matched, onquery, onreference, hardclip, softclip, notcigar,
};


read_type get_rg_type_from_str(const std::string &str);
cigar_character_type what_is_this_cigar(char c);

#define CIGAR_CHARACTERS "DHIMNPSX="
#define TAB "\t"



class sam_tag{
public:
    std::string tag;
    char type;
    std::string value;
    sam_tag() : tag(""),type('X'),value(""){}
    sam_tag(const std::string &tag, char type, const std::string &value) : tag(tag),type(type),value(value){
    }
    sam_tag(const std::string &whole_string) : tag(whole_string.substr(0,2)),type(whole_string[3]),value(whole_string.substr(5)){
    }
};

class sam_record{
public:

    //SAM fields
    std::string  qname;
    int     flag;
    std::string  rname;
    int     pos;
    int  mapq;
    std::string  cig_str; 
    cigar   cgr;
    std::string  mate_rname;
    int     mate_pos;
    int     tlen;
    std::string  seq;
    std::string  qual; 
    std::map<std::string, sam_tag> tags;

    //Utility Fields
    int head_clip_range;
    int tail_clip_range;
    bool is_fully_mapped; 
    int end_pos;
    read_type rg;

    int insertion_pos;

    bool invalid = false;
    
    sam_record( const std::string &line);

    void analyze_cigar();
    void set_rg( read_type rg);
    friend std::ostream& operator<<(std::ostream& os, const sam_record& sr);
};


std::ostream& operator<<(std::ostream& os, const sam_record& sr);

char *get_nrv();
void reverseComplete (char *seq, char *rcSeq , int length);

#endif
