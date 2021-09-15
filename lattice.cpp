#include "lattice.h"
#include "vectorMath.h"
#include <iostream>
#include <math.h>

const double pi = 3.14159265259;

Lattice::Lattice(std::vector<double> & v1, std::vector<double> & v2, std::vector<double> & v3, double a)
{
    x = v1;
    y = v2;
    z = v3;
    alpha = a;
    calcVolume();
    updateRecip();
    calcMadelung();
}

void Lattice::updateRecip()
{
    a = vecMath::scalarMult(2.0*pi/volume,vecMath::cross(y,z));
    b = vecMath::scalarMult(2.0*pi/volume,vecMath::cross(z,x));
    c = vecMath::scalarMult(2.0*pi/volume,vecMath::cross(x,y));
}

void Lattice::calcVolume()
{
    std::vector<double> yz = vecMath::cross(y,z);
    volume = abs(vecMath::dot(x,yz));
}

double Lattice::realSpace_Madelung(int nmax)
{
    double val = 0.0;
    for (int nx = -nmax; nx <= nmax; nx++)
    {
        for (int ny = -nmax; ny <= nmax; ny++)
        {
            for (int nz = -nmax; nz <= nmax; nz++)
            {
                if ((nx==0)&&(ny==0)&&(nz==0))
                {
                    continue;
                }
                std::vector<double> n = vecMath::sum(vecMath::scalarMult(nx,x),vecMath::scalarMult(ny,y),vecMath::scalarMult(nz,z));
                double norm = vecMath::norm(n);
                val += erfc(alpha*norm)/norm;
            }
        }
    }
    return val;
}

double Lattice::recipSpace_Madelung(int nmax)
{
    double val = 0.0;
    for (int na = -nmax; na <= nmax; na++)
    {
        for (int nb = -nmax; nb <= nmax; nb++)
        {
            for (int nc = -nmax; nc <= nmax; nc++)
            {
                if ((na==0)&&(nb==0)&&(nc==0))
                {
                    continue;
                }
                std::vector<double> n = vecMath::sum(vecMath::scalarMult(na,a),vecMath::scalarMult(nb,b),vecMath::scalarMult(nc,c));
                double norm = vecMath::norm(n);
                val += exp(-norm*norm/(4*alpha*alpha))/(norm*norm);
            }
        }
    }
    return (4.0*pi/volume)*val;

}

void Lattice::calcMadelung()
{
    //Real Space
    double prev = 0;
    double real = 0;
    for (int n = 1; n < 100; n++)
    {
        real = realSpace_Madelung(n);
        if (abs(real-prev)< 1.0e-6)
        {
            real_nmax = n;
            break;
        }
        prev = real;
    }

    //Recip Space
    prev = 0;
    double recip = 0.0;
    for (int n = 1; n < 100; n++)
    {
        recip = recipSpace_Madelung(n);
        if (abs(recip-prev)<1.0e-6)
        {
            recip_nmax = n;
            break;
        }
        prev = recip;
    }
    double wo_term = real + recip - 2*alpha/sqrt(pi);
    double w_term = wo_term - pi/(alpha*alpha*volume);
    std::cout << alpha << " " << wo_term << " " << w_term << std::endl;
    madelung = real + recip - 2*alpha/sqrt(pi) - pi/(alpha*alpha*volume);
}

double Lattice::ewald(const std::vector<double> & r)
{
    double real = 0.0;
    //real space
    for (int nx=-real_nmax; nx<=real_nmax; nx++)
    {
        for (int ny=-real_nmax; ny<=real_nmax; ny++)
        {
            for (int nz=-real_nmax; nz<=real_nmax; nz++)
            {
                std::vector<double> n = vecMath::sum(vecMath::scalarMult(nx,x),vecMath::scalarMult(ny,y),vecMath::scalarMult(nz,z));
                std::vector<double> rn = vecMath::sum(r,n);
                double norm = vecMath::norm(rn);
                real += erfc(alpha*norm)/norm;
            }
        }
    }
    //recip space
    double recip = 0.0;
    for (int na=-recip_nmax; na<=recip_nmax; na++)
    {
        for (int nb=-recip_nmax; nb<=recip_nmax; nb++)
        {
            for (int nc=-recip_nmax; nc<=recip_nmax; nc++)
            {
                if ((na==0)&&(nb==0)&&(nc==0))
                {
                    continue;
                }
                std::vector<double> g = vecMath::sum(vecMath::scalarMult(na,a),vecMath::scalarMult(nb,b),vecMath::scalarMult(nc,c));
                double norm = vecMath::norm(g);
                double dot = vecMath::dot(g,r);
                recip += exp(-norm*norm/(4*alpha*alpha))*cos(dot)/(norm*norm);
            }
        }
    }
    recip *= 4*pi/volume;

    return real + recip; //- pi/(volume*alpha*alpha);

}

double Lattice::C3d()
{
    double prefactor = 0.25 * pow(volume,4./3.);
    double c3d;

    double alpha = 1;
    while (alpha > 1e-5)
    {
        double recip = 0.0;
        double prev = 1e10;
        for (int nmax = 1; nmax <= 50; nmax++)
        {
            double val = 0.0;
            for (int na = -nmax; na <= nmax; na++)
            {
                for (int nb = -nmax; nb <= nmax; nb++)
                {
                    for (int nc = -nmax; nc <= nmax; nc++)
                    {
                        if ((na==0)&&(nb==0)&&(nc==0))
                        {
                            continue;
                        }
                        std::vector<double> n = vecMath::sum(vecMath::scalarMult(na,a),vecMath::scalarMult(nb,b),vecMath::scalarMult(nc,c));
                        double norm = vecMath::norm(n);
                        val += 4*pi/volume * norm * exp(-alpha * norm * norm);
                    }
                }
            }
            if (abs(val-prev)<1.0e-8)
            {
                recip = val;
                break;
            }
            prev = val;
            if (nmax == 50)
                std::cerr << " error" << std::endl;
        }
    
        c3d = prefactor * (1.0/(pi*alpha*alpha) - recip);
        std::cout << alpha << " " << c3d << std::endl;
        alpha *= 0.9;
    }

    return c3d;
}
