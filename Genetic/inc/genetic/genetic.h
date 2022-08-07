#ifndef INC_GENETIC_GENETIC
#define INC_GENETIC_GENETIC

#include <algorithm>
#include <cassert>
#include <numeric>
#include <span>
#include <concepts>
#include <vector>
#include <random>
#include <iostream>

#include <genetic/opt.h>
#include <genetic/operations.h>
#include <genetic/tournament.h>

namespace g
{
    struct Solution
    {
        std::vector<double> x;
        double value;
    };

    template<typename F>
    Solution solve(F & f, const std::span<double> x0,
                   const std::span<double> low, const std::span<double> high,
                   const Options & opt, const int iters)
    {
        const int nGenes = x0.size();
        assert(nGenes == low.size());
        assert(nGenes == high.size());

        std::vector<double> scales(nGenes);
        for (int i = 0; i < nGenes; i++)
        {
            scales[i] = (high[i] - low[i]) / 2.0; 
        }

        std::vector<double> tmp(nGenes);
        auto fNorm = [&](const std::span<double> x)
        {
            for (int i = 0; i < nGenes; i++)
            {
                tmp[i] = scales[i] * (x[i] + 1.0) + low[i];
            }
            return f(tmp);
        };

        std::vector<double> x0Norm(nGenes);
        for (int i = 0; i < nGenes; i++)
        {
            x0Norm[i] = (x0[i] - low[i]) / scales[i] - 1.0;
        }

        auto sol = solveNormalized(fNorm, x0Norm, opt, iters);
        for (int i = 0; i < nGenes; i++)
        {
            sol.x[i] = scales[i] * (sol.x[i] + 1.0) + low[i];
        }
        return sol;
    }

    template<typename F>
    Solution solveNormalized(F & f, const std::span<double> x0, 
                             const Options & opt, const int iters)
    {
        const int nGenes = x0.size();
        const int n = opt.populationSize;

        std::vector<std::vector<double>> p0(n);
        for (auto & g : p0) g.resize(nGenes);
        std::vector<std::vector<double>> p1(n);
        std::vector<int> idx(n);
        std::iota(idx.begin(), idx.end(), 0);

        std::default_random_engine rng(std::random_device{}());
        std::uniform_int_distribution<int> idxDist(0, n - 1);
        std::uniform_int_distribution<int> crossoverDist(0, nGenes - 2);
        std::uniform_real_distribution<double> pDist(0, 1);
        std::normal_distribution<double> smallN(0, opt.small.sigma);
        std::normal_distribution<double> bigN(0, opt.big.sigma);

        // Generate initial population
        std::copy_n(x0.begin(), nGenes, p0[0].begin());
        for (int i = 1; i < n; i++) 
        {
            uniformMutation(p0[i], rng);
        }

        auto * ptrParents = &p0;
        auto * ptrChildren = &p1;
        std::vector<double> values(n);
        for (int i = 0; i < iters; i++)
        {
            auto & parents = *ptrParents;
            auto & children = *ptrChildren;

            // Selection
            std::transform(parents.begin(), parents.end(), values.begin(), [&](auto & v) { return f(v); });
            for (int j = 0; j < n; j++)
            {
                const int winner = runTournament(values, opt.tournamentSize, rng, idxDist);
                children[j] = parents[winner];
            }

            // Mutation
            for (int j = 0; j < n; j++)
            {
                auto & curr = children[j];
                if (pDist(rng) <= opt.pUniformMutation)
                {
                    uniformMutation(curr, rng);
                }
                // Big + small mutations
                for (int k = 0; k < nGenes; k++)
                {
                    if (pDist(rng) <= opt.big.p)
                    {
                        const double delta = bigN(rng);
                        curr[k] += delta;
                    }
                    if (pDist(rng) <= opt.small.p)
                    {
                        const double delta = smallN(rng);
                        curr[k] += delta;
                    }
                    curr[k] = std::clamp(curr[k], -1.0, 1.0);
                }
            }

            // Crossovers
            std::shuffle(idx.begin(), idx.end(), rng);
            for (int j = 0; j < n - 1; j += 2)
            {
                if (pDist(rng) <= opt.pCrossover)
                {
                    const int cIdx = crossoverDist(rng);
                    simpleCrossover(children[idx[j]], children[idx[j + 1]], cIdx);
                }
            }

            // Find best parent and copy to output
            auto bestIt = std::min_element(values.begin(), values.end());
            children[0] = parents[bestIt - values.begin()];

            std::cout << i << "\t" << *bestIt << "\n";

            std::swap(ptrParents, ptrChildren);
        }

        // Evaluate one last time
        auto & final = *ptrParents;
        std::transform(final.begin(), final.end(), values.begin(), [&](auto & v) { return f(v); });
        auto bestIt = std::min_element(values.begin(), values.end());

        Solution sol;
        sol.value = *bestIt;
        sol.x = final[bestIt - values.begin()];

        return sol; 
    }
}

#endif
