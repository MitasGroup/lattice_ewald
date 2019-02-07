#ifndef VECTORMATH_H
#define VECTORMATH_H

#include <vector>

namespace vecMath {
    double dot(const std::vector<double> & a, const std::vector<double> & b);
    std::vector<double> cross(const std::vector<double> & a, const std::vector<double> & b);
    double norm(const std::vector<double> & a);
    std::vector<double> sum(const std::vector<double> & a,const std::vector<double> & b);
    std::vector<double> sum(const std::vector<double> & a,const std::vector<double> & b,const std::vector<double> & c);
    template<class T> 
    std::vector<double> scalarMult(T a, const std::vector<double> & v)
    {
        std::vector<double> mult = v;
        for (int d = 0; d < mult.size(); d++)
            mult[d] *= a;
        return mult;
    }
    template<class T> 
    std::vector<double> scalarMult(const std::vector<double> & v,T a)
    {
        return vecMath::scalarMult(a,v);
    }
};


#endif
