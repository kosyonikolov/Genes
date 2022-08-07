#include <iostream>
#include <span>
#include <cmath>

#include <genetic/opt.h>
#include <genetic/genetic.h>

double f1(const std::span<double> v)
{
    const double a = 10;

    const double x = v[0];
    const double y = v[1];

    double val = 2 * a + x * x + y * y - a * std::cos(2 * M_PI * x) - a * std::cos(2 * M_PI * y);
    return val;
}

int main()
{
    std::array<double, 2> low = {-10, -10};
    std::array<double, 2> high = { 10, 10};
    std::array<double, 2> x0 = {5, 5};

    g::Options opt;
    opt.populationSize = 15;
    opt.tournamentSize = 3;
    opt.pCrossover = 0.5;
    opt.pUniformMutation = 0.02;
    opt.big.p = 0.2;
    opt.big.sigma = 0.2;
    opt.small.p = 0.7;
    opt.small.sigma = 0.05;

    auto sol = g::solve(f1, x0, low, high, opt, 100);

    std::cout << sol.value << "\t" << sol.x[0] << ", " << sol.x[1] << "\n";
}
