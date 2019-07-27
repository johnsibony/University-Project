/*
 Ce fichier permet de regrouper en une même classe les différentes méthode de pricing (mêmes méthodes que pour le swap de variance) d'une option européenne et d'afficher leur prix
 */


#ifndef europeanOptionPricer_hpp
#define europeanOptionPricer_hpp

#include "montecarloPrice.hpp"
#include <string>

class EuropeanOptionPricer
{
    
public:
    
    EuropeanOptionPricer();
    EuropeanOptionPricer(double T, int N, double sO, double vO, double K); // hyperparameters constructor
    EuropeanOptionPricer(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO, double K, double threshold); // all parameters constructor
    double function_price(int nb, double k);
    double get_analytic_price(int nbSim);
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
    double _K;
    double _threshold;
    EulerScheme _euler;
    TGScheme _TG;
    QEScheme _QE;
    
};

#endif
