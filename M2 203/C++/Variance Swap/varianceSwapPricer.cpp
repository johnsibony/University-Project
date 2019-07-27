#include "utilsClass.hpp"
#include "varianceSwapPricer.hpp"
#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;



VarianceSwapPricer::VarianceSwapPricer(): _kappa(0.5), _rho(-0.9), _sigma(1), _theta(0.04), _r(0.0), _T(10), _N(500), _sO(100), _vO(0.04), _threshold(1.5), _analytic(_T, _N, _vO, _kappa, _rho, _sigma, _theta, _r), _euler(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO), _TG(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO), _QE(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO, _threshold)
{
    
}
        
VarianceSwapPricer::VarianceSwapPricer(double T, int N, double sO, double vO): _T(T), _N(N), _sO(sO),_vO(vO), _analytic(_T, _N, 0.04, 1.5, -0.5, 2, 0.02, 0.05), _euler(_T, _N, 1.5, 2, 0.02, 0.05, -0.5, _sO, _vO), _TG(_T, _N, 1.5, 2, 0.02, 0.05, -0.5, _sO, _vO), _QE(_T, _N, 1.5, 2, 0.02, 0.05, -0.5, _sO, _vO, 1.5)
{
    
}
    
VarianceSwapPricer::VarianceSwapPricer(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO, double threshold): _T(T), _N(N), _kappa(kappa), _sigma(sigma), _theta(theta), _r(r), _rho(rho), _sO(sO), _vO(vO), _threshold(threshold), _analytic(_T, _N, _vO, _kappa, _rho, _sigma, _theta, _r), _euler(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO), _TG(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO), _QE(_T, _N, _kappa, _sigma, _theta, _r, _rho, _sO, _vO, _threshold)
{
        
}
    
My_complex VarianceSwapPricer::get_analytic_price()
{
    My_complex price (0, 0);
        
    for (int i=2; i<=_N; i++)
    {
        price += _analytic.variance_asset(i);
    }
    price = (1 / _T) * (_analytic.variance_asset_initial() + price) * (-(pow(100, 2)));
        
    return price;
}
    
double VarianceSwapPricer::get_analytic_price(string complexe_part)
{
    My_complex price (0, 0);
        
    for (int i=2; i<=_N; i++)
    {
        price += _analytic.variance_asset(i);
    }
        
    price = (1 / _T) * (_analytic.variance_asset_initial() + price) * (-(pow(100, 2)));
        
    if(complexe_part == "re")
    {
        return price.real();
    }
    else if(complexe_part == "im")
    {
        return price.imag();
    }
    else
    {
        throw string("Complexe part name unrecognized");
    }
}
    
double VarianceSwapPricer::get_montecarlo_price(string scheme, int nbMC)
{
    double priceMonteCarlo = 0;
    int pass = 0;  // number of MC iteration returning nan
    double tps = clock();
    
    ofstream pricefile("/Users/johnsibony/desktop/swapprice"+ scheme +".txt");
    ofstream timefile("/Users/johnsibony/desktop/swaptime"+ scheme +".txt");
        
    for (int j=1; j<=nbMC; j++)
    {
        double price = 0;
        double * logAssetsPrice = 0;
            
        if(scheme == "euler")
        {
            logAssetsPrice = _euler.log_asset_simulation();
        }
        else if(scheme == "tg")
        {
            logAssetsPrice = _TG.log_asset_simulation();
        }
        else if(scheme == "qe")
        {
            logAssetsPrice = _QE.log_asset_simulation();
        }
        else
        {
            throw string("Scheme name unrecognized");
        }
            
        for (int i=1; i<=_N; i++)
        {
            price += pow(logAssetsPrice[i] - logAssetsPrice[i-1], 2);
        }

        if( isnan((pow(100, 2) / _T) * price) == 1)  // if the pricing does not fit, we pass (value explosion due to approximiation)
        {
            pass ++;
            continue;
        }
        
        priceMonteCarlo += (pow(100, 2) / _T) * price;
    
        pricefile << priceMonteCarlo / (j-pass) << endl;
        timefile << (clock() - tps) / (double)CLOCKS_PER_SEC << endl;
        //cout << "MonteCarlo price :" << priceMonteCarlo / (j-pass) << endl;
        //cout << "time " << (clock() - tps) / (double)CLOCKS_PER_SEC << endl;
    }
        
    return priceMonteCarlo / (nbMC - pass);
}
