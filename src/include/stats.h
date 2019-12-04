#ifndef __PAMIR_STATS_H__
#define __PAMIR_STATS_H__

#include <utility>
#include <cmath>
#include <numeric>
#include <iterator>
#include <functional>

template<class ITER>
static inline std::pair<double, double> get_distribution(ITER first, ITER last) {
    double mean = 0;
    size_t count = std::distance(first, last);
    long long sum = 0;
    sum = std::accumulate(first, last, sum);
    mean = static_cast<double>(sum) / count;

    double variance = 0;
    for (auto iter = first; iter != last; ++iter) {
        variance = variance + ((*iter - mean) * (*iter - mean));
    }
    variance = variance / count;

    return std::make_pair(mean, std::sqrt(variance));
}

template<class ITER, typename ACCESS>
static inline std::pair<double, double> get_distribution(ITER first, ITER last, ACCESS accessor_func) {

    size_t count = std::distance(first, last);
    long long sum = 0;
    for (auto iter = first; iter != last; ++iter) {
        sum += accessor_func(iter);
    }
    double mean = static_cast<double>(sum) / count;

    double variance = 0;
    for (auto iter = first; iter != last; ++iter) {
        variance = variance + ((accessor_func(iter) - mean) * (accessor_func(iter) - mean));
    }
    variance = variance / count;

    return std::make_pair(mean, std::sqrt(variance));
}

#endif