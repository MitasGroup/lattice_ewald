#ifndef LATTICE_H
#define LATTICE_H

#include <vector>

class Lattice 
{
    public:
        Lattice(std::vector<double> & v1, std::vector<double> & v2, std::vector<double> & v3,double a);
        double volume;
        double madelung;
        std::vector<double> x,y,z;
        std::vector<double> a,b,c;
        double ewald(const std::vector<double> & r);
        double C3d();
    private:
        int real_nmax;
        int recip_nmax;
        double alpha;
        void updateRecip();
        void calcVolume();
        void calcMadelung();
        double realSpace_Madelung(int nmax);
        double recipSpace_Madelung(int nmax);
};

#endif
