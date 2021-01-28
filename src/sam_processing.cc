#include "sam_processing.h"
#include "cigar.h"

#include <cstring>
#include <iostream>
#include <sstream>

using std::stringstream;
using std::map;
using std::pair;
using std::string;
using std::map;
using std::vector;


std::array<std::string, 7> read_types{{
    "BreakPointSupport",
    "BreakPointSupportOEA",
    "BreakPointSupportOrphan",
    "ReferenceSupport",
    "FalseSupport",
    "NoSupport",
    "Unset",
}};

string &get_read_type_str(read_type type) {
    return read_types[static_cast<int>(type)];
}

read_type get_rg_type_from_str(const string &str) {
    for (int i = 0; i < read_types.size(); i++) {
        if (str == read_types[i]) {
            return static_cast<read_type>(i);
        }
    }
    return read_type::unset;
}

cigar_character_type what_is_this_cigar(char c) {
    switch (c) {
        case 'M':
        case '=':
        case 'X':
            return cigar_character_type::matched;
        case 'D':
        case 'N':
            return cigar_character_type::onreference;
        case 'I':
        case 'P':
            return cigar_character_type::onquery;
        case 'H':
            return cigar_character_type::hardclip;
        case 'S':
            return cigar_character_type::softclip;
        default:
            return cigar_character_type::notcigar;
    }
}


cigar::cigar(const cigar &cig) : cigarray(cig.cigarray) {
}

cigar::cigar(const vector <pair<int, char>> &cigarray) : cigarray(cigarray) {
}

cigar::cigar(const string &cgr) {
    size_t prev = -1;
    size_t index = cgr.find_first_of(CIGAR_CHARACTERS);
    while (index != string::npos) {
        int len = std::stoi(cgr.substr(1 + prev, index));
        cigarray.push_back(std::make_pair(len, cgr[index]));
        prev = index;
        index = cgr.find_first_of(CIGAR_CHARACTERS, index + 1);
    }
}


std::ostream &operator<<(std::ostream &os, const cigar &cig) {
    if (cig.cigarray.size() == 0) {
        os << "*";
        return os;
    }
    for (auto iter = cig.cbegin(); iter != cig.cend(); ++iter) {
        os << iter->first << iter->second;
    }
    return os;
}

decltype(cigar::cigarray.cbegin()) cigar::cbegin() const {
    return cigarray.cbegin();
}

decltype(cigar::cigarray.begin()) cigar::begin() {
    return cigarray.begin();
}


decltype(cigar::cigarray.cend()) cigar::cend() const {
    return cigarray.cend();
}

decltype(cigar::cigarray.rbegin()) cigar::rbegin() {
    return cigarray.rbegin();
}

decltype(cigar::cigarray.rend()) cigar::rend() {
    return cigarray.rend();
}

decltype(cigar::cigarray.end()) cigar::end() {
    return cigarray.end();
}

auto cigar::operator[](size_t index) const {
    return cigarray[index];
}


void sam_record::analyze_cigar() {
    end_pos = pos;
    head_clip_range = 0;
    tail_clip_range = 0;
    insertion_pos = pos;
    bool is_first = true;
    is_fully_mapped = true;
    if (cig_str == "*")
        is_fully_mapped = false;
    for (auto cc : cgr) {  // cc is a pair of <int, char>
        auto type = what_is_this_cigar(cc.second);
        if (type == cigar_character_type::matched) {
            end_pos += cc.first;
        } else if (type == cigar_character_type::onquery) {
            insertion_pos = end_pos;
            is_fully_mapped = false;
        } else if (type == cigar_character_type::onreference) {
            end_pos += cc.first;
            is_fully_mapped = false;
        } else if (type == cigar_character_type::softclip) {
            is_fully_mapped = false;
            if (is_first) {
                head_clip_range = cc.first;
            } else {
                tail_clip_range = cc.first;
            }
        }
        is_first = false;
    }
}

std::ostream &operator<<(std::ostream &os, const sam_record &sr) {
    os << sr.qname << TAB
       << sr.flag << TAB
       << sr.rname << TAB
       << sr.pos << TAB
       << sr.mapq << TAB
       << sr.cgr << TAB
       << sr.mate_rname << TAB
       << sr.mate_pos << TAB
       << sr.tlen << TAB
       << sr.seq << TAB
       << sr.qual;
    for (auto pair : sr.tags) {
        sam_tag &tag = pair.second;
        os << TAB << tag.tag << ":" << tag.type << ":" << tag.value;
    }
    return os;
}

void sam_record::set_rg(read_type rg) {
    this->rg = rg;
    sam_tag rg_tag("RG:Z:" + read_types[static_cast<int>(rg)]);
    tags[rg_tag.tag] = rg_tag;
}

sam_record::sam_record(const string &line) {
    stringstream sp(line);

    sp >> qname
       >> flag
       >> rname
       >> pos
       >> mapq
       >> cig_str
       >> mate_rname
       >> mate_pos
       >> tlen
       >> seq
       >> qual;
    string next;
    while (sp >> next) {
        sam_tag _tag(next);
        tags[_tag.tag] = _tag;
    }
    auto rg_index = tags.find("RG");
    if (rg_index != tags.end()) {
        rg = get_rg_type_from_str(rg_index->second.value);
    } else {
        rg = read_type::nosup;
    }
    cgr = cigar(cig_str);
    analyze_cigar();
}
/**
 * Generates the reverse complement lookup table
 * @return Reverse complement lookup table
 */
char *get_nrv() {
    static char *nRev = NULL;
    if (nRev == NULL) {
        nRev = (char *) malloc(128);
        memset(nRev, 'N', 128);
        nRev['A'] = 'T';
        nRev['C'] = 'G';
        nRev['T'] = 'A';
        nRev['G'] = 'C';
    }
    return nRev;
}

/**
 * Reverse complement the sequence (from mrsfast source code)
 * @param seq Original sequence
 * @param rcSeq Reverse complement sequence
 * @param length Length of the sequence
 */
void reverseComplete(char *seq, char *rcSeq, int length) {
    int i;
    char *nRev = get_nrv();
    seq += length - 1;
    for (i = 0; i < length; i++) {
        rcSeq[i] = nRev[*(seq--)];
    }
    rcSeq[length] = '\0';
}
