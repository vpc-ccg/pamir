#ifndef CIGAR_CIGAR
#define CIGAR_CIGAR
#include <string>
#include <iostream>
#include <vector>
#include <limits>

class cigar{
    std::vector<std::pair<int,char>> cigarray;
    public:



    cigar() {}
    cigar(const std::string &cgr); 
    cigar(const cigar &cgr); 
    cigar(const std::vector<std::pair<int,char>> &cgr_vec); 

    decltype(cigarray.cbegin()) cbegin() const;
    decltype(cigarray.cend()) cend() const;

    decltype(cigarray.begin()) begin();
    decltype(cigarray.rbegin()) rbegin();
    decltype(cigarray.rend()) rend();
    decltype(cigarray.end()) end();
    template <class iter>
    auto refill( const iter beg, const iter end){
        cigarray.clear();
        return cigarray.insert(cigarray.end(),beg,end);
    }
    auto operator [](size_t index) const;
    auto back() const;
    friend std::ostream& operator<<(std::ostream& os, const cigar& cig);
};

#endif
