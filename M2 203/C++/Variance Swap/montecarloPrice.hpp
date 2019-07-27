/*
 Ce fichier permet de simuler les log prix de l'asset S selon la méthode de Monte Carlo pour les schémas d'Euler, TG et QE
 */

#ifndef montecarloPrice_hpp
#define montecarloPrice_hpp

#include "utilsClass.hpp"


/*
 Class to compute the next stochastic asset price with its stochastic variance
*/

class HestonModel
{
    
    public:
    
    HestonModel(double T, int N, double kappa, double rho, double sigma, double theta, double r);
    double S(double gauss, double s, double v_current, double v_next); // asset price simulation given the current one with Broade-Kaya scheme
    double* chi_moments(double v); // compute m^2 and s^2
    double time_path(); // we must define the time path to compute the stochastic asset

    protected: // The inheritated class of HestonModel could use its private atributes
    
    double _T;
    int _N;
    double _kappa;
    double _rho;
    double _sigma;
    double _r;
    double _theta;
    RandomTool random;
    
};

/*
 Class to compute the next stochastic asset price using the Euler scheme
*/

class EulerScheme: public HestonModel
{
    
    public:
    
    EulerScheme(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO);
    double* log_asset_simulation();
    double S(double gauss, double s, double v); // redefinition of the asset price simulation (Euler scheme of instead of Broade-Kaya scheme for S)
    double V(double gauss, double v); // stochastic variance given the current one
    
    private:
    
    double _sO;
    double _vO;
    
};

/*
 Class to compute the next stochastic asset price using the TG scheme
*/

class TGScheme: public HestonModel
{
    
    public:
    
    TGScheme(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO);
    double* log_asset_simulation();
    double V(double gauss, double v);
    double* TG_moments(double v);
    double newton_search(double phi);
    double equation(double r, double phi);
    double equation_derive(double r, double phi);
    double density(double x);
    double density_derive(double x);
    double CDF(double x);
    
    private:
    
    double _sO;
    double _vO;
    
};

/*
 Class to compute the next stochastic asset price using the QE scheme
*/

class QEScheme: public HestonModel
{
    
    public:
    
    QEScheme(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO, double threshold);
    double* log_asset_simulation();
    double V(double gauss, double v);
    double next_below(double m, double Z, double phi);
    double* compute_a_b(double m, double phi);
    double next_above(double m, double U, double phi);
    double* compute_p_beta(double m, double phi);
    double inverse_cumulative(double U, double p, double beta);
    
    private:
    
    double _sO;
    double _vO;
    double _threshold;
    
};

#endif
