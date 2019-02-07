#include "vectorMath.h"
#include <assert.h>
#include <math.h>

double vecMath::dot(const std::vector<double> & a, const std::vector<double> & b)
{
    assert(a.size() == b.size());
    double val = 0.0;
    for (int d = 0; d < a.size(); d++)
    {
        val += a[d]*b[d];
    }
    return val;
}

std::vector<double> vecMath::cross(const std::vector<double> & a, const std::vector<double> & b)
{
    assert(a.size() == b.size());
    assert(a.size() == 3);
    std::vector<double> cross(3);
    for (int d = 0; d < 3; d++)
    {
        cross[d] = a[(d+1)%3]*b[(d+2)%3] - a[(d+2)%3]*b[(d+1)%3];
    }
    return cross;
}

double vecMath::norm(const std::vector<double> & a)
{
    return sqrt(vecMath::dot(a,a));
}

std::vector<double> vecMath::sum(const std::vector<double> & a, const std::vector<double> & b)
{
    assert(a.size() == b.size());
    std::vector<double> res(a.size());
    for (int d = 0; d < a.size(); d++)
    {
        res[d] = a[d]+b[d];
    }
    return res;
}

std::vector<double> vecMath::sum(const std::vector<double> & a, const std::vector<double> & b, const std::vector<double> & c)
{
    assert(a.size() == b.size());
    assert(b.size() == c.size());
    std::vector <double> ab = vecMath::sum(a,b);
    return vecMath::sum(ab,c);
}

