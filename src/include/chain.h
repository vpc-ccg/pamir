/*
 * clasp is written by Christian Otto, Bioinformatics, University of Leipzig
 * https://www.bioinf.uni-leipzig.de/Software/clasp/
 *
 * This is a wrapper class designed for Pamir.
 */

#ifndef CHAIN_H
#define CHAIN_H

#include <vector>

#include "clasp/manopt.h"

using namespace std;

#define SOP     ((unsigned char) (0 << 0))
#define LIN     ((unsigned char) (1 << 0))

struct MaxChainInfo {
    pair<int, int> qrange;
    pair<int, int> rrange;
    float score;
    int len;
};

struct seed {
    int qPos;
    int sPos;
    int len;
};

class ClaspChain {
    private:
        char chainmode;
        double lambda;
        double epsilon;
        int maxgap;

    public:
        ClaspChain(char chainmode = SOP, double lambda = 0.5, double epsilon = 0, int maxgap = -1);
        MaxChainInfo get_max_chain(vector<seed> seeds);
};

#endif