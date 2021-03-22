#include "chain.h"

#include <chrono>

#ifdef __cplusplus
extern "C" {
#endif

#include "clasp/container.h"
#include "clasp/sltypes.h"
#include "clasp/slchain.h"

#ifdef __cplusplus
}
#endif

#include <iostream>

//#include "clasp/sltypes.h"
//#include "clasp/slchain.h"
//#include "clasp/container.h"

ClaspChain::ClaspChain(char chainmode, double lambda, double epsilon, int maxgap) :
                    chainmode(chainmode), lambda(lambda), epsilon(epsilon), maxgap(maxgap) {}

MaxChainInfo ClaspChain::get_max_chain(vector<seed>& seeds) {
    Container *fragments;

    int num;

    MaxChainInfo max_chain;

    max_chain.score = -1;

    /* initialization */
    fragments = (Container *) malloc(sizeof(Container));
    bl_containerInit(fragments, 1000, sizeof(slmatch_t));

    /* put the fragments in the container */
    slmatch_t frag;
    for (int i = 0; i < seeds.size(); i++) {
        bl_slmatchInit(&frag, 0);
        frag.p = seeds[i].sPos;
        frag.i = seeds[i].qPos;
        frag.q = frag.j = seeds[i].len;
        frag.scr = seeds[i].len;
        bl_containerAdd(fragments, &frag);
    }

    /* sort fragments */
    auto t1 = chrono::high_resolution_clock::now();
    qsort(fragments->contspace, bl_containerSize(fragments),
          sizeof(slmatch_t), cmp_slmatch_qsort);
    auto t2 = chrono::high_resolution_clock::now();
    sort_time += std::chrono::duration<double, std::milli>(t2-t1).count();

    int begin = 0;
//    for (int i = 1; i <= bl_containerSize(fragments); i++) {
        /*
        * end of fragments list or different database sequence
        * --> process fragment[begin]...fragment[i-1], write output
        *     and free chains (less memory consumption with large input files)
        */
//        if (i == bl_containerSize(fragments) ||
//            ((slmatch_t *) bl_containerGet(fragments, begin))->subject !=
//            ((slmatch_t *) bl_containerGet(fragments, i))->subject) {

//            if (chainmode == SOP) {
                /* only use chaining without clustering if no ids are specified */
                //bl_slChainSop((slmatch_t *) info.fragments->contspace + begin, i - begin,
                //        info.epsilon, info.lambda);
                auto t11 = chrono::high_resolution_clock::now();
//                bl_slClusterSop((slmatch_t *) fragments->contspace + begin, i - begin, epsilon, lambda, maxgap);
                bl_slClusterSop((slmatch_t *) fragments->contspace, seeds.size(), epsilon, lambda, maxgap);
                auto t22 = chrono::high_resolution_clock::now();
                chain_time += std::chrono::duration<double, std::milli>(t22-t11).count();
//            }
//            else {
//                //bl_slChainLin((slmatch_t *) info.fragments->contspace + begin, i - begin,
//                //        info.epsilon, info.lambda);
//                bl_slClusterLin((slmatch_t *) fragments->contspace + begin, i - begin, epsilon, lambda, maxgap);
//            }

//            for (int j = begin; j < i; j++) {
            for (int j = 0; j < seeds.size(); j++) {
                slmatch_t *match = (slmatch_t *) bl_containerGet(fragments, j);

                if (match->chain) {
                    slchain_t *chain = (slchain_t *) match->chain;

                    if (chain->scr > max_chain.score) {
//                        std::cerr << chain->scr << std::endl;
                        max_chain.score = chain->scr;
                        max_chain.len = chain->q;
                        max_chain.qrange = {chain->i, chain->i + chain->j - 1};
                        max_chain.rrange = {chain->p, chain->p + chain->q - 1};
                        max_chain.gaps_size = abs((max_chain.qrange.second - max_chain.qrange.first) -
                                (max_chain.rrange.second - max_chain.rrange.first));
                    }

                    bl_slchainDestruct(chain);
                    free(chain);
                    match->chain = NULL;
                }
            }
//            begin = i;
//        }
//    }

    if (fragments) {
        for (int i = 0; i < bl_containerSize(fragments); i++) {
            slmatch_t *sl = (slmatch_t *) bl_containerGet(fragments, i);
            if (sl->chain != NULL) {
                exit(-1);
                bl_slchainDestruct(sl->chain);
                free(sl->chain);
            }
        }
        bl_containerDestruct(fragments, bl_slmatchDestruct);
        free(fragments);
    }

    return max_chain;
}