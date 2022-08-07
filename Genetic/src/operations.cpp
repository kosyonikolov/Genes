#include <cassert>
#include <genetic/operations.h>

namespace g
{
    void simpleCrossover(std::span<double> a, std::span<double> b, const int idx)
    {
        const int n = a.size();
        assert(b.size() == n);
        assert(idx >= 0 && idx < n);

        if (idx < n / 2)
        {
            for (int i = 0; i <= idx; i++) std::swap(a[i], b[i]);
        }
        else
        {
            for (int i = idx + 1; i < n; i++) std::swap(a[i], b[i]);
        }
    }
}