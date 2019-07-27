/*
 Ce fichier permet de regrouper en une même classe les différentes méthode de pricing d'un swap de variance et d'afficher leur prix
 */


#ifndef varianceSwapPricer_hpp
#define varianceSwapPricer_hpp

#include "analyticPrice.hpp"
#include "montecarloPrice.hpp"
#include <string>



/*
 Class to compute the Variance swap price according different scheme
*/

class VarianceSwapPricer
{
    
    public:
    
    VarianceSwapPricer();
    VarianceSwapPricer(double T, int N, double sO, double vO); // hyperparameters constructor
    VarianceSwapPricer(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO, double threshold); // all parameters constructor
    My_complex get_analytic_price(); // return the complex analytic price
    double get_analytic_price(std::string complexe_part); // return the real or imaginary part of the complex analytic price
    double get_montecarlo_price(std::string scheme, int nbMC); // return the monte-carlo price given the scheme type and the number of simulation

    private:
    
    double _T;
    int _N;
    double _kappa;
    double _sigma; // analytic notation ('epsilon' name for MonteCarlo assumed to be named 'sigma')
    double _r;
    double _theta;
    double _sO;
    double _vO;
    double _rho;
    double _threshold;
    Analytic _analytic;
    EulerScheme _euler;
    TGScheme _TG;
    QEScheme _QE;
    
};


#endif
