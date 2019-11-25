#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <sstream>
#include <array>
#include <algorithm>

#include <cstdlib>
#include <cstring>

#include "sam_processing.h"
#include "cigar.h"

#include <edlib.h>

using std::string;
using std::vector;
using std::array;
using std::stringstream;
using std::map;
using std::pair;


namespace cram_process {

    enum class sup_type {
        rightbp, leftbp, both, ref, nosup,
    };

    sup_type pad_ref_support(sam_record &rec, int flank_len, int insert_len) {

        auto iter = rec.cgr.begin();
        int pos = rec.pos;
        int cigar_len;
        int cigar_char;
        vector <pair<int, char>> new_cigar_array;
        for (; iter != rec.cgr.end(); ++iter) {
            cigar_len = iter->first;
            cigar_char = iter->second;
            auto type = what_is_this_cigar(cigar_char);
            if (type == cigar_character_type::onreference || type == cigar_character_type::matched) {
                if (pos + cigar_len > flank_len) {
                    break;
                }
                pos += cigar_len;
            }
            new_cigar_array.push_back(std::make_pair(cigar_len, cigar_char));
        }

        if (flank_len - pos > 0) {
            new_cigar_array.push_back(std::make_pair(flank_len - pos, cigar_char));
        }
        ++iter;
        while (iter != rec.cgr.end() && what_is_this_cigar(iter->second) == cigar_character_type::onreference) {
            insert_len += iter->first;
            ++iter;
        }
        new_cigar_array.push_back(std::make_pair(insert_len, 'N'));

        if (pos + cigar_len - flank_len > 0) {
            new_cigar_array.push_back(std::make_pair(pos + cigar_len - flank_len, cigar_char));
        }

        for (; iter != rec.cgr.end(); ++iter) {
            cigar_len = iter->first;
            cigar_char = iter->second;
            new_cigar_array.push_back(std::make_pair(cigar_len, cigar_char));
        }
        rec.cgr = new_cigar_array;

        return sup_type::ref;
    }

#define INSERT_MIN_CONSIDER_PERCENT 0.6

    bool cigar_contains_insertion(const cigar &cgr, int insert_len) {

        double min_insert_len = INSERT_MIN_CONSIDER_PERCENT * insert_len;
        int insert_count = 0;
        for (auto iter = cgr.cbegin(); iter != cgr.cend(); ++iter) {
            auto pair = *iter;
            auto type = what_is_this_cigar(pair.second);
            if (type == cigar_character_type::onquery) {
                insert_count += pair.first;
            }
        }
        return insert_count >= min_insert_len;
    }

    sup_type remap_if_on_breakpoint(sam_record &rec, char *seq, int flank_len, int insert_len) {
        if (rec.pos < flank_len + 1 && rec.end_pos > flank_len + 1) {
            if (rec.head_clip_range == 0 && rec.tail_clip_range == 0) {
                if (!cigar_contains_insertion(rec.cgr, insert_len)) {
                    rec.set_rg(read_type::refsup);
                    return pad_ref_support(rec, flank_len, insert_len);

                }
            }
        }
        if (rec.pos - rec.head_clip_range > flank_len + insert_len || rec.end_pos + rec.tail_clip_range < flank_len) {
            return sup_type::nosup;
        }

        rec.set_rg(read_type::bpsup);

        int ref_len = 2 * flank_len + insert_len;
        EdlibAlignResult result = edlibAlign(rec.seq.c_str(), rec.seq.size(), seq, ref_len,
                                             edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

        int new_pos = result.startLocations[0];


        string new_cigar_str = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
        cigar new_cgr(new_cigar_str);

        rec.cgr = new_cgr;
        rec.cig_str = new_cigar_str;
        rec.pos = 1 + new_pos;
        rec.end_pos = 1 + result.endLocations[0];
        if (rec.pos < flank_len && rec.end_pos > flank_len + insert_len) {
            return sup_type::both;
        } else if (rec.pos < flank_len) {
            return sup_type::rightbp;
        } else if (rec.end_pos > insert_len + flank_len) {
            return sup_type::leftbp;
        }
        return sup_type::nosup;
    }

    void recalibrate(sam_record &rec, int ins_loci, int flank_len, string insert_id) {
        rec.rname = insert_id;
        rec.pos = rec.pos - ins_loci + flank_len;
        rec.end_pos = rec.end_pos - ins_loci + flank_len;
    }

    void shift_if_right_flank(sam_record &rec, int ins_loci, int flank_len, int ins_len) {
        if (rec.pos < ins_loci) {
            return;
        }
        rec.pos += ins_len;
        rec.end_pos += ins_len;
    }

    void cut_right_boundary(sam_record &rec, int ins_loci, int flank_len) {
        if (rec.end_pos < ins_loci + flank_len) {
            return;
        }
        int move_len = rec.end_pos - ins_loci - flank_len;

        auto cigar_iter = rec.cgr.rbegin();

        vector <pair<int, char>> new_cigar_array;
        int cigar_len;
        int cigar_char;

        if (rec.tail_clip_range > 0) {
            ++cigar_iter;
            if (cigar_iter == rec.cgr.rend()) {
                return;
            }
        }
        new_cigar_array.push_back(std::make_pair(rec.tail_clip_range + move_len, 'S'));


        while (move_len != 0) {
            cigar_len = cigar_iter->first;
            cigar_char = cigar_iter->second;
            auto type = what_is_this_cigar(cigar_char);
            if (type == cigar_character_type::matched) {

                int min_len = std::min(cigar_len, move_len);
                cigar_len -= move_len;
                rec.end_pos -= min_len;
                if (cigar_len <= 0) {
                    move_len = -cigar_len;
                    ++cigar_iter;
                    continue;
                } else {
                    new_cigar_array.push_back(std::make_pair(cigar_len, cigar_char));
                    ++cigar_iter;
                    break;
                }
            } else if (type == cigar_character_type::onreference) {
                if (what_is_this_cigar(new_cigar_array.back().second) == cigar_character_type::softclip) {
//                    move_len -= cigar_len;
                    //    rec.end_pos -= cigar_len;
                    if (move_len <= 0) {
                        ++cigar_iter;
                        break;
                    }
                } else {
                    std::cerr << "Check up : " << rec.cig_str << std::endl;
                }
            } else if (type == cigar_character_type::onquery) {
                if (what_is_this_cigar(new_cigar_array.back().second) == cigar_character_type::softclip) {
                    pair<int, char> last_pair = new_cigar_array.back();
                    new_cigar_array.pop_back();
                    new_cigar_array.push_back(std::make_pair(last_pair.first + cigar_len, 'S'));
                } else {
                    std::cerr << "Check up : " << rec.cig_str << std::endl;
                }
            } else if (type == cigar_character_type::softclip) {

            }
            if (cigar_iter == rec.cgr.rend()) {
                if (what_is_this_cigar(new_cigar_array.back().second) == cigar_character_type::softclip &&
                    move_len > 0) {
                    pair<int, char> last_pair = new_cigar_array.back();
                    new_cigar_array.pop_back();
                    new_cigar_array.push_back(std::make_pair(last_pair.first + move_len, 'S'));
                }
                break;
            }

            ++cigar_iter;
        }

        for (; cigar_iter != rec.cgr.rend(); ++cigar_iter) {
            new_cigar_array.push_back(*cigar_iter);
        }

        rec.cgr.refill(new_cigar_array.rbegin(), new_cigar_array.rend());

    }

    void make_unmapped(sam_record &rec) {
        rec.pos = -1;
        rec.cig_str = "*";
        rec.rname = "*";
        rec.flag &= 4073;
        rec.flag |= 4;
        rec.invalid = true;
    }

    void cut_left_boundary(sam_record &rec, int ins_loci, int flank_len) {
        if (rec.pos > ins_loci - flank_len) {
            return;
        }
        int move_len = 1 + ins_loci - flank_len - rec.pos;
        auto cigar_iter = rec.cgr.begin();

        vector <pair<int, char>> new_cigar_array;
        int cigar_len;
        char cigar_char;

        if (rec.head_clip_range > 0) {
            ++cigar_iter;
            if (cigar_iter == rec.cgr.end()) {
                return;
            }
        }
        if (move_len + rec.tail_clip_range + rec.head_clip_range >= rec.seq.size()) {
            make_unmapped(rec);
            return;
        }

        new_cigar_array.push_back(std::make_pair(rec.head_clip_range + move_len, 'S'));

        while (move_len != 0) {
            cigar_len = cigar_iter->first;
            cigar_char = cigar_iter->second;
            auto type = what_is_this_cigar(cigar_char);
            if (type == cigar_character_type::matched) {
                int min_len = std::min(cigar_len, move_len);
                cigar_len -= move_len;

                rec.pos += min_len;
                if (cigar_len <= 0) {
                    move_len = -cigar_len;
                    ++cigar_iter;
                    continue;
                } else {
                    new_cigar_array.push_back(std::make_pair(cigar_len, cigar_char));
                    ++cigar_iter;
                    break;
                }

            } else if (type == cigar_character_type::onreference) {
                if (what_is_this_cigar(new_cigar_array.back().second) == cigar_character_type::softclip) {
//    move_len -= cigar_len;
                    rec.pos += cigar_len;
                    if (move_len <= 0) {
                        if (cigar_iter != rec.cgr.end()) {
                            ++cigar_iter;
                        }
                        break;
                    }
                } else {
                    std::cerr << "Check up : " << rec.cig_str << std::endl;
                }
            } else if (type == cigar_character_type::onquery) {
                if (what_is_this_cigar(new_cigar_array.back().second) == cigar_character_type::softclip) {

                    pair<int, char> last_pair = new_cigar_array.back();
                    new_cigar_array.pop_back();
                    new_cigar_array.push_back(std::make_pair(last_pair.first + cigar_len, 'S'));
                } else {
                    std::cerr << "Check up : " << rec.cig_str << std::endl;
                }
            } else if (type == cigar_character_type::softclip) {

            }

            if (cigar_iter == rec.cgr.end()) {
                break;
            }

            ++cigar_iter;
        }


        for (; cigar_iter != rec.cgr.end(); ++cigar_iter) {
            new_cigar_array.push_back(*cigar_iter);
        }

        rec.cgr = new_cigar_array;

    }


    sup_type remap_unmapped(sam_record &read, const string &chr_name, char *seq, int flank_len, int insert_len) {
        int seq_len = 2 * flank_len + insert_len;
        int edlib_K = 10;
        int reversed = 0;
        EdlibAlignResult result = edlibAlign(read.seq.c_str(), read.seq.size(), seq, seq_len,
                                             edlibNewAlignConfig(edlib_K, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
        if (result.numLocations == 0) {
            char *revseq = (char *) malloc(read.seq.size());
            reverseComplete((char *) read.seq.c_str(), revseq, read.seq.size());
            result = edlibAlign(revseq, read.seq.size(), seq, seq_len,
                                edlibNewAlignConfig(edlib_K, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

            if (result.numLocations == 0) {
                return sup_type::nosup;
            }
            read.seq = revseq;
            reversed = 16;
        }

        read.pos = 1 + result.startLocations[0];
        read.end_pos = 1 + result.endLocations[0];
        string new_cigar_str = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
        read.cig_str = new_cigar_str;
        read.cgr = cigar(new_cigar_str);

        read.rname = chr_name;
        read.mapq = 33;

        unsigned mask = static_cast<unsigned>(-1);
        mask ^= 4;
        read.flag &= mask;

        read.flag |= reversed;
        if (read.pos < flank_len && read.end_pos > flank_len + insert_len) {
            read.set_rg(read_type::bpsup);
            return sup_type::both;
        }
        if (read.pos < flank_len && read.end_pos > flank_len) {
            read.set_rg(read_type::bpsup);
            return sup_type::leftbp;
        }
        if (read.pos > flank_len && read.end_pos > flank_len + insert_len) {
            read.set_rg(read_type::bpsup);
            return sup_type::rightbp;
        }
        if (read.pos > flank_len && read.end_pos < flank_len + insert_len) {
            read.set_rg(read_type::bpsupoea);
            return sup_type::nosup;
        }
        read.set_rg(read_type::nosup);
        return sup_type::nosup;
    }

    read_type read_type_more_important(read_type a, read_type b) {

        if (a == read_type::bpsup || b == read_type::bpsup) {
            return read_type::bpsup;
        }
        if (a == read_type::bpsupoea || b == read_type::bpsupoea) {
            return read_type::bpsupoea;
        }
        if (a == read_type::refsup || b == read_type::refsup) {
            return read_type::refsup;
        }
        if (a == read_type::falsesup || b == read_type::falsesup) {
            return read_type::falsesup;
        }
        return read_type::nosup;
    }

    void increment_supports(sup_type sup, double &left, double &right, double &ref) {
        switch (sup) {
            case sup_type::both:
                left += 1;
                right += 1;
                break;
            case sup_type::leftbp:
                left += 1;
                break;
            case sup_type::rightbp:
                right += 1;
                break;
            case sup_type::ref:
                ref += 1;
                break;
            default:
                break;
        }
    }

    int main(int argc, char **argv) {

        int insert_loci = atoi(argv[1]);
        int insert_len = atoi(argv[2]);
        int flank_len = atoi(argv[3]);
        int min_frag_len = atoi(argv[4]);
        int max_frag_len = atoi(argv[5]);
        string insert_id = argv[6];
        char *fasta = argv[7];
        int range_begin = insert_loci - flank_len;
        int range_end = insert_loci + flank_len;

        string prev_line = "-1";
        string prev_qname = "-1";
        map <string, vector<sam_record>> reads;

        double left_sup = 0, right_sup = 0, ref_sup = 0;

        for (std::string line; std::getline(std::cin, line);) {
            sam_record read(line);
            if (read.rg == read_type::unset) {
                read.set_rg(read_type::nosup);
            }
            cut_left_boundary(read, insert_loci, flank_len);
            cut_right_boundary(read, insert_loci, flank_len);
            shift_if_right_flank(read, insert_loci, flank_len, insert_len);
            recalibrate(read, insert_loci, flank_len, insert_id);

            if (read.invalid) {
                continue;
            }

            if ((read.flag & 4) == 4) {
                read.flag &= 69;
                remap_unmapped(read, insert_id, fasta, flank_len, insert_len);
            } else {
                sup_type support = remap_if_on_breakpoint(read, fasta, flank_len, insert_len);
                increment_supports(support, left_sup, right_sup, ref_sup);
            }

            reads[read.qname].push_back(read);
        }
        std::cout << left_sup << "\t" << right_sup << "\t" << ref_sup << "\n";
        for (auto p : reads) { // p is a pair<string, vector<string>>
            vector <sam_record> &read_pair = p.second;
            if (read_pair.size() == 1) { //Single End
                sam_record &read = read_pair[0];
                read.flag = (read.flag & 20);
                read.mate_rname = "*";
                read.mate_pos = 0;
                read.tlen = 0;
                std::cout << read << "\n";
            } else {

                read_pair[0].mate_pos = read_pair[1].pos;
                read_pair[1].mate_pos = read_pair[0].pos;

                read_pair[0].mate_rname = read_pair[1].rname;
                read_pair[1].mate_rname = read_pair[0].rname;

                read_type new_rg = read_type_more_important(read_pair[0].rg, read_pair[1].rg);
                read_pair[0].set_rg(new_rg);
                read_pair[1].set_rg(new_rg);


                int min_pos = std::min({read_pair[0].pos, read_pair[0].end_pos, read_pair[1].pos, read_pair[1].end_pos},
                                       [](int a, int b) { return a < b; });
                int max_pos = std::max({read_pair[0].pos, read_pair[0].end_pos, read_pair[1].pos, read_pair[1].end_pos},
                                       [](int a, int b) { return a < b; });

                int new_tlen = max_pos - min_pos;

                read_pair[0].tlen = new_tlen;
                read_pair[1].tlen = new_tlen;

                unsigned mask = static_cast<unsigned>(-1);
                mask ^= 8;
                if ((read_pair[0].flag & 4) == 0) {
                    read_pair[1].flag &= mask;
                }
                if ((read_pair[1].flag & 4) == 0) {
                    read_pair[0].flag &= mask;
                }
                std::cout << read_pair[0] << "\n";
                std::cout << read_pair[1] << "\n";


            }
        }
        return 0;
    }
}
