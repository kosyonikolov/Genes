#ifndef INC_GENETIC_OPT
#define INC_GENETIC_OPT

namespace g
{
    struct MutationOptions
    {
        double sigma;
        double p; // per gene probability
    };

    struct Options
    {
        int populationSize;
        int tournamentSize;
        MutationOptions small, big;
        double pUniformMutation;
        double pCrossover;
    };
}

#endif
