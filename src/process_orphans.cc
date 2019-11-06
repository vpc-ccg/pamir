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


namespace orphan_process{
#define INSERT_MIN_CONSIDER_PERCENT 0.6
    bool cigar_contains_insertion(const cigar &cgr, int insert_len){

        double min_insert_len = INSERT_MIN_CONSIDER_PERCENT * insert_len;
        int insert_count = 0;
        for(auto iter = cgr.cbegin(); iter!=cgr.cend();++iter){
            auto pair = *iter;
            auto type = what_is_this_cigar(pair.second);
            if( type == cigar_character_type::onquery){
                insert_count += pair.first;
            }
        }
        return insert_count >= min_insert_len;
    }


    int main(int argc, char **argv){

        string prev_line = "-1";
        string prev_qname  = "-1";
        map<string, vector<sam_record>> reads;
        for (std::string line; std::getline(std::cin, line);) {
            if(line[0] == '@'){ //Header
                std::cout << line << "\n";
                continue;
            }
            sam_record read(line);
            read.set_rg(read_type::bpsuporphan);

            if((read.flag & 4)  == 4  || (read.flag &8)  == 8){ //We skip if mate unmapped
                continue;
            }
             
            reads[read.qname].push_back(read);
        }

        for(auto p : reads){ // p is a pair<string, vector<string>>

            vector<sam_record> &read_pair = p.second;
            
            if(read_pair.size() == 1){
                std::cerr << "Something is wrong, no pair for : " << read_pair[0].qname << std::endl;
            }else{
                bool flag = true;
                for( auto r  : read_pair){ // skip if there are any soft clips
                    if ( r.head_clip_range !=0){
                        flag = false;
                    }
                    if ( r.tail_clip_range !=0){
                        flag = false;
                    }
                    if (cigar_contains_insertion(r.cgr,24)){
                            flag = false;        
                    }
                }
                if ( read_pair[0].rname != read_pair[1].rname){
                    flag = false;
                }
                if( flag){
                    std::cout << read_pair[0] << "\n";
                    std::cout << read_pair[1] << "\n";
                }

            }
        }
        return 0;
    }
}
