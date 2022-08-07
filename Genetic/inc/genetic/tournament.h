#ifndef INC_GENETIC_TOURNAMENT
#define INC_GENETIC_TOURNAMENT

#include <cassert>
#include <span>

namespace g
{
    template<typename RNG, typename Distribution>
    int runTournament(const std::span<double> v, const int size, RNG & rng, Distribution & dist)
    {
        assert(size > 0);
        
        int best = dist(rng);
        double bestVal = v[best];
        for (int i = 1; i < size; i++)
        {
            int cand = dist(rng);
            double candVal = v[cand];
            if (candVal < bestVal)
            {
                best = cand;
                bestVal = candVal;
            }
        }

        return best;
    }
}

#endif
