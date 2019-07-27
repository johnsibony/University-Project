#include "europeanOptionPricer.hpp"
#include "utilsClass.hpp"
#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;


EuropeanOptionPricer::EuropeanOptionPricer(): _kappa(1.5), _rho(-0.5), _sigma(1), _theta(0.02), _r(0.0), _T(1), _N(100), _sO(100), _vO(0.04), _K(70), _threshold(1.5), _euler(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO), _TG(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO), _QE(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO, _threshold)
{
    
}

EuropeanOptionPricer::EuropeanOptionPricer(double T, int N, double sO, double vO, double K): _T(T), _N(N), _sO(sO),_vO(vO), _K(K), _euler(_T, _N, 1.5, 2, 0.02, 0.05, -0.5, _sO, _vO), _TG(_T, _N, 1.5, 2, 0.02, 0.05, -0.5, _sO, _vO), _QE(_T, _N, 1.5, 2, 0.02, 0.05, -0.5, _sO, _vO, 1.5)
{
    
}

EuropeanOptionPricer::EuropeanOptionPricer(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO, double K, double threshold): _T(T), _N(N), _kappa(kappa), _sigma(sigma), _theta(theta), _r(r), _rho(rho), _sO(sO), _vO(vO), _K(K), _threshold(threshold), _euler(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO), _TG(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO), _QE(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO, _threshold)
{
    
}

double EuropeanOptionPricer::function_price(int nb, double k)
{
    My_complex j(0.0,1.0);
    My_complex alpha(0.0,0.0);
    My_complex beta(0.0,0.0);
    if(nb == 1)
    {
        beta = _kappa - (_rho * _sigma) - (_rho * _sigma * k * j);
        alpha = -pow(k,2)- 0.5*k*j + j*k;
    }
    if(nb == 2)
    {
        beta = _kappa - 0*(_rho * _sigma) - (_rho * _sigma * k * j);
        alpha = -pow(k,2)- 0.5*k*j + 0*j*k;
    }
    
    My_complex d = sqrt(pow(beta,2) - 2*alpha*pow(_sigma,2));
    My_complex r_ = (beta - d) / pow(_sigma,2);
    My_complex rpos = (beta + d) / pow(_sigma,2);
    My_complex g = r_ / rpos;
    My_complex D = ((1 - exp(-d * _T)) * r_) / (pow(_sigma, 2) * (1 - (g * exp(-d * _T))));
    My_complex C = _kappa * (r_ * _T - ((2/pow(_sigma,2)) * log((1 - g * exp(-d * _T))/(1 - g))));
    
    return (exp(k * log((_sO*exp(_r*_T))/_K) * j + (C * _theta) + (D * _vO)) / (k * j)).real();
}

double EuropeanOptionPricer::get_analytic_price(int nbSim)
{
    double borne_inf = 0;
    double borne_sup = 100;
    double step = (borne_sup-borne_inf)/nbSim;
    double riemann_integral1 = 0;
    double riemann_integral2 = 0;
    
    for (int i = 1; i <= nbSim; i ++)
    {
        riemann_integral1 += function_price(1, borne_inf + i*step) * step; // sum up each small rectangle
        riemann_integral2 += function_price(2, borne_inf + i*step) * step;
    }
    return _sO * (0.5 + ((1/M_PI) * riemann_integral1)) - _K * exp(-_r * _T) * (0.5 + ((1/M_PI) * riemann_integral2));
}

double EuropeanOptionPricer::get_montecarlo_price(string scheme, int nbMC)
{
    double price = 0;
    double * logAssetsPrice = 0;
    double tps = clock();
    
    ofstream pricefile("/Users/johnsibony/desktop/optionprice"+ scheme +".txt");
    ofstream timefile("/Users/johnsibony/desktop/optiontime"+ scheme +".txt");
        
    if(scheme == "euler")
    {
        for (int i = 1; i <= nbMC; i ++)
        {
        logAssetsPrice = _euler.log_asset_simulation();
        price += max(exp(logAssetsPrice[_N]) - _K, 0.0);
            
        pricefile << price / i << endl;
        timefile << (clock() - tps) / (double)CLOCKS_PER_SEC << endl;
        //cout << price / i << endl;
        //cout << (clock() - tps) / (double)CLOCKS_PER_SEC << endl;
        }
        return price / nbMC;
    }
    else if(scheme == "tg")
    {
        for (int i = 1; i <= nbMC; i ++)
        {
            logAssetsPrice = _TG.log_asset_simulation();
            price += max(exp(logAssetsPrice[_N]) - _K, 0.0);
            
            pricefile << price / i << endl;
            timefile << (clock() - tps) / (double)CLOCKS_PER_SEC << endl;
            //cout << price / i << endl;
            //cout << (clock() - tps) / (double)CLOCKS_PER_SEC << endl;
        }
        return price / nbMC;
    }
    else if(scheme == "qe")
    {
        for (int i = 1; i <= nbMC; i ++)
        {
            logAssetsPrice = _QE.log_asset_simulation();
            price += max(exp(logAssetsPrice[_N]) - _K, 0.0);
            
            pricefile << price / i << endl;
            timefile << (clock() - tps) / (double)CLOCKS_PER_SEC << endl;
            //cout << price / i << endl;
            //cout << (clock() - tps) / (double)CLOCKS_PER_SEC << endl;
        }
        return price / nbMC;
    }
    else
    {
        throw string("Scheme name unrecognized");
    }
    
}
