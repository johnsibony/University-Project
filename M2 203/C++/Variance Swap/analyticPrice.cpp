#include "analyticPrice.hpp"
#include "utilsClass.hpp"

using namespace std;



All_functions::All_functions(double kappa, double rho, double sigma): _kappa(kappa), _rho(rho), _sigma(sigma)
{
        
}

My_complex All_functions::a(double omega)
{
    My_complex j(0,1);
    return _kappa - _rho * _sigma * omega * j;
}
    
My_complex All_functions::b(double omega)
{
    My_complex j(0,1);
    return sqrt(pow(a(omega), 2) + pow(_sigma, 2) * (pow(omega, 2) + omega * j));
}
    
My_complex All_functions::g(double omega)
{
    return (a(omega) + b(omega)) / (a(omega) - b(omega));
}





C::C(double kappa, double rho, double sigma, double theta, double r) : _theta(theta), _r(r), All_functions(kappa, rho, sigma)
{
        
}
    
My_complex C::value_C(double tau, double omega)
{
    My_complex j(0,1);
    
    if(isnan(g(omega).imag()) == 1) // when a = b, g is not defined (g = inf), so we rewrite C te deal with thus limit.
    {
        return _r * tau * (- 1 + omega * j)  + ((_kappa * _theta) / (pow(_sigma, 2))) * (tau * (a(omega) + b(omega)) - 2 * tau * b(omega));
    }
    else
    {
        return _r * tau * (- 1 + omega * j)  + ((_kappa * _theta) / (pow(_sigma, 2))) * (tau * (a(omega) + b(omega)) - 2 * log((1 - g(omega) * exp(tau * b(omega))) / (1-g(omega))));
    }
}
    
My_complex C::first_derivative(double tau, double h)
{
    return ( value_C(tau, h) - value_C(tau, 0) ) / h;
}
    
My_complex C::second_derivative(double tau, double h) // compute the second derivative of C with the Taylor Young formula
{
    return ( value_C(tau, h) - 2 * value_C(tau, 0) + value_C(tau, -h) ) / pow(h, 2);
}





D::D(double kappa, double rho, double sigma): All_functions(kappa, rho, sigma)
{
        
}
    
My_complex D::value_D(double tau, double omega) // no problem compared with C : even if g = inf, the normal formula of D can be used and D wille be zero
{
    My_complex j(0,1);
    return ((a(omega) + b(omega)) / pow(_sigma, 2)) * (1 - exp(-tau * b(omega))) / (1 - g(omega) * exp(-tau * b(omega)));
}
    
My_complex D::first_derivative(double tau, double h)
{
    return ( value_D(tau, h) - value_D(tau, 0) ) / h;
}
    
My_complex D::second_derivative(double tau, double h) // compute the second derivative of D with the Taylor Young formula
{
    return ( value_D(tau, h) - 2 * value_D(tau, 0) + value_D(tau, -h) ) / pow(h, 2);
}





Analytic::Analytic(double T, int N, double vO, double kappa, double rho, double sigma, double theta, double r): _kappa(kappa), _rho(rho), _sigma(sigma), _r(r), _theta(theta), C(kappa, rho, sigma, theta, r), D(kappa, rho, sigma), _T(T), _N(N), _vO(vO)
{
        
}
    
double Analytic::time(int i)
{
    return i * _T / _N;
}
    
double Analytic::time_path()
{
    return double(_T)/_N;
}
    
double Analytic::c(int i)
{
    return (2 * _kappa) / (pow(_sigma, 2) * (1 - exp(-_kappa * time(i-1))));
}
    
double Analytic::W(int i)
{
    return c(i) * _vO * exp(- _kappa * time(i-1));
}
    
double Analytic::q()
{
    return (2 * _kappa * _theta) / pow(_sigma, 2);
}
    
My_complex Analytic::variance_asset_initial()
{
    double h = 0.000001;
    return pow(D.first_derivative(time_path(), h), 2) * pow(_vO, 2) + ( 2 * C.first_derivative(time_path(), h) * D.first_derivative(time_path(), h) - D.second_derivative(time_path(), h) ) * _vO + pow(C.first_derivative(time_path(), h), 2) - C.second_derivative(time_path(), h);
}
    
My_complex Analytic::variance_asset(int i)
{
    double h = 0.000001;
    return pow(D.first_derivative(time_path(), h), 2) * ( q() + 2 * W(i) + pow(q() + W(i), 2) / pow(c(i), 2) ) + ( 2 * C.first_derivative(time_path(), h) * D.first_derivative(time_path(), h) - D.second_derivative(time_path(), h) ) * ( (q() + W(i)) / c(i) ) + pow(C.first_derivative(time_path(), h), 2) - C.second_derivative(time_path(), h);
}
