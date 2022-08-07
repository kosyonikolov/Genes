#ifndef INC_GENETIC_OPERATIONS
#define INC_GENETIC_OPERATIONS

#include <algorithm>
#include <span>
#include <random>

namespace g
{
    template<typename RNG>
    void uniformMutation(std::span<double> v, RNG & rng)
    {
        std::uniform_real_distribution<double> dist(-1, 1);
        std::generate(v.begin(), v.end(), [&](){ return dist(rng); });
    }

    void simpleCrossover(std::span<double> a, std::span<double> b, const int idx);
}

#endif
