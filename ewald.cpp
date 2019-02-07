#include <iostream>
#include <vector>
#include "lattice.h"
#include <iomanip>
#include "vectorMath.h"
#include <fstream>
#include <string>

int main(int argc, char ** argv) 
{
    std::cout.precision(10);

    std::vector<double> x(3);
    std::vector<double> y(3);
    std::vector<double> z(3);

    if (argc < 2)
    {
        std::cout << "Usage: ./ewald --latvec x1 x2 x3 y1 y2 y3 z1 z2 z3 [options]" << std::endl;
        std::cout << "options: " << std::endl;
        std::cout << "  --alpha a " << std::endl;
        std::cout << "  --plot_ewald" << std::endl;
        std::cout << "  --madelung_scaling" << std::endl;
        return 0;
    }

    int iargc = 1;
    double alpha = 1.0;
    bool plot_ewald=false;
    bool madelung_scaling=false;
    while(iargc < argc)
    {
        std::string a(argv[iargc]);
        if (a == "--latvec")
        {
            x[0] = atof(argv[++iargc]);
            x[1] = atof(argv[++iargc]);
            x[2] = atof(argv[++iargc]);
            y[0] = atof(argv[++iargc]);
            y[1] = atof(argv[++iargc]);
            y[2] = atof(argv[++iargc]);
            z[0] = atof(argv[++iargc]);
            z[1] = atof(argv[++iargc]);
            z[2] = atof(argv[++iargc]);
        }
        else if (a == "--alpha")
        {
            alpha = atof(argv[++iargc]);
        }
        else if (a == "--plot_ewald")
        {
            plot_ewald = true;
        }
        else if (a == "--madelung_scaling")
        {
            madelung_scaling = true;
        }
        else 
        {
            std::cout << "Error: Invalid argument" << std::endl;
            return 0;
        }
        iargc++;
    }

    Lattice lat(x,y,z,alpha);
    std::cout << "Cell Volume: " << lat.volume << std::endl;
    std::cout << "Madelung   : " << lat.madelung << std::endl;

    if (madelung_scaling)
    {
        std::ofstream mad;
        mad.open("madelung.dat");
        for (double scale = 1.0; scale <= 10; scale *= 1.1)
        {
            std::vector<double> sx = vecMath::scalarMult(scale,lat.x);
            std::vector<double> sy = vecMath::scalarMult(scale,lat.y);
            std::vector<double> sz = vecMath::scalarMult(scale,lat.z);
            Lattice slat(sx,sy,sz,alpha);
            std::cout << "Scale: " << scale << std::endl;
            mad << slat.volume << " " << slat.madelung << std::endl;
        }
        mad.close();

    }

    if (plot_ewald)
    {
        std::ofstream fv1,fv2,fv3,fv4,fv5,fv6,fv7;
        fv1.open("coulomb_ewald_100.dat");
        fv2.open("coulomb_ewald_010.dat");
        fv3.open("coulomb_ewald_001.dat");
        fv4.open("coulomb_ewald_110.dat");
        fv5.open("coulomb_ewald_101.dat");
        fv6.open("coulomb_ewald_011.dat");
        fv7.open("coulomb_ewald_111.dat");
        for (double pv = 0.01; pv <= 1.0; pv+=0.005)
        {
            std::vector<double> v1 = vecMath::scalarMult(pv,lat.x);
            std::vector<double> v2 = vecMath::scalarMult(pv,lat.y);
            std::vector<double> v3 = vecMath::scalarMult(pv,lat.z);
            std::vector<double> v4 = vecMath::scalarMult(pv,vecMath::sum(lat.x,lat.y));
            std::vector<double> v5 = vecMath::scalarMult(pv,vecMath::sum(lat.x,lat.z));
            std::vector<double> v6 = vecMath::scalarMult(pv,vecMath::sum(lat.y,lat.z));
            std::vector<double> v7 = vecMath::scalarMult(pv,vecMath::sum(lat.x,lat.y,lat.z));
            fv1 << vecMath::norm(v1) << " " << lat.ewald(v1) << std::endl;
            fv2 << vecMath::norm(v2) << " " << lat.ewald(v2) << std::endl;
            fv3 << vecMath::norm(v3) << " " << lat.ewald(v3) << std::endl;
            fv4 << vecMath::norm(v4) << " " << lat.ewald(v4) << std::endl;
            fv5 << vecMath::norm(v5) << " " << lat.ewald(v5) << std::endl;
            fv6 << vecMath::norm(v6) << " " << lat.ewald(v6) << std::endl;
            fv7 << vecMath::norm(v7) << " " << lat.ewald(v7) << std::endl;
        }
        fv1.close();
        fv2.close();
        fv3.close();
        fv4.close();
        fv5.close();
        fv6.close();
        fv7.close();
    }

}
