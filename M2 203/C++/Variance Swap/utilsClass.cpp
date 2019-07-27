#include "utilsClass.hpp"

using namespace std;



My_complex::My_complex(double a, double b): complex<double>(a, b)
{
        
}
    
My_complex::My_complex(complex<double> a): complex<double>(a)
{
        
}

My_complex operator* (const double nb, const My_complex z)
{
    return My_complex(nb * z.real(), nb * z.imag());
}

My_complex operator+ (const double nb,  const My_complex z)
{
    return My_complex(nb + z.real(), z.imag());
}

My_complex operator- (const double nb,  const My_complex z)
{
    return My_complex(nb - z.real(), -z.imag());
}

My_complex operator* (const int nb,  const My_complex z)
{
    return My_complex(nb * z.real(), nb * z.imag());
}

My_complex operator+ (const int nb,  const My_complex z)
{
    return My_complex(nb + z.real(), z.imag());
}

My_complex operator- (const int nb,  const My_complex z)
{
    return My_complex(nb - z.real(), -z.imag());
}





RandomTool::RandomTool(): _rho(0) 
{
        
}
    
RandomTool::RandomTool(double rho): _rho(rho)
{
    if(rho > 1 || rho < -1)  // We raise an error if the correlation coefficient does not belong to [-1, 1]
    {
        throw string("The correlation coefficients must belong to [-1, 1]");
    }
}
    
double RandomTool::randomUniform()
{
    return  rand()/(double)RAND_MAX;
}
    
double* RandomTool::boxMuller()  // We return a pointer to an array 
{
    double U1 = randomUniform();
    double U2 = randomUniform();
    const double u = sqrt(-2 * log(U1)) * cos(2 * M_PI * U2);
    const double v = sqrt(-2 * log(U1)) * sin(2 * M_PI * U2);
    
    double *tab = new double [2];
    tab[0] = u;
    tab[1] = (_rho * u) + sqrt(1 - pow(_rho, 2)) * v;
        
    return tab;
}
