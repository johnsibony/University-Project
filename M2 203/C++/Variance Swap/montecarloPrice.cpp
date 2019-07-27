#include "montecarloPrice.hpp"
#include <iostream>

using namespace std;



HestonModel::HestonModel(double T, int N, double kappa, double rho, double sigma, double theta, double r): random(rho), _T(T), _N(N), _kappa(kappa), _sigma(sigma), _theta(theta), _r(r)
{
    
}

double HestonModel::S(double gauss, double s, double v_current, double v_next) // asset next price given the current s one
{
    double K0 = -(_rho * _kappa * _theta * time_path()) / _sigma;
    double K1 = 0.5 * time_path() * ((_kappa * _rho / _sigma) - 0.5) -_rho/_sigma;
    double K2 = 0.5 * time_path() * ((_kappa * _rho / _sigma) - 0.5) +_rho/_sigma;
    double K3 = 0.5 * time_path() * (1 - pow(_rho, 2));
    double K4 = 0.5 * time_path() * (1 - pow(_rho, 2));
    
    return s + K0 + K1*v_current + K2*v_next + sqrt(K3*v_current + K4*v_next) * gauss;
}

double* HestonModel::chi_moments(double v) // compute m^2 and s^2
{
    double *moment = new double [2];
    
    moment[0] = pow(_theta + (v - _theta) * exp(-_kappa * time_path()), 2);
    moment[1] = ( (v * pow(_sigma, 2) * exp(-_kappa * time_path())) / (_kappa) ) * (1 - exp(-_kappa * time_path())) + (_theta * pow(_sigma, 2)/(2 * _kappa)) * pow(1 - exp(-_kappa * time_path()), 2);
    
    return moment;
}

double HestonModel::time_path()  // constant path. Maybe the lecturer would change it.
{
    return double(_T)/_N;
}





EulerScheme::EulerScheme(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO): HestonModel(T, N, kappa, rho, sigma, theta, r), _sO(sO),  _vO(vO)
{
    
}

double* EulerScheme::log_asset_simulation()
{
    double * simulationAsset = new double [_N+1];
    double * simulationVariance = new double [_N+1];
    simulationAsset[0] = log(_sO);  // BE CAREFULL : we take the log initial price
    simulationVariance[0] = _vO;
    
    for (int i=1; i<=_N; i++)
    {
        double *a = random.boxMuller(); // a couple of correlated gaussian law
        simulationAsset[i] = S(a[0], simulationAsset[i-1], simulationVariance[i-1]);
        simulationVariance[i] = V(a[1], simulationVariance[i-1]);
    }
    
    return simulationAsset;
}

double EulerScheme::S(double gauss, double s, double v)
{
    return s - (0.5 * max(v, 0.0) * time_path()) + (sqrt(max(v, 0.0)) * gauss * sqrt(time_path()));
}

double EulerScheme::V(double gauss, double v)
{
    return v + (_kappa * (_theta - max(v, 0.0)) * time_path() ) + _sigma * sqrt(max(v, 0.0)) * gauss * sqrt(time_path());
}





TGScheme::TGScheme(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO): HestonModel(T, N, kappa, rho, sigma, theta, r), _sO(sO),  _vO(vO)
{
    
}

double* TGScheme::log_asset_simulation()
{
    double gauss;
    double * simulationAsset = new double [_N+1];
    double * simulationVariance = new double [_N+1];
    simulationAsset[0] = log(_sO); // BE CAREFULL : we take the log initial price
    simulationVariance[0] = _vO;
    
    for (int i=1; i<=_N; i++)
    {
        gauss = random.boxMuller()[0]; //first gaussian law
        simulationVariance[i] = V(gauss, simulationVariance[i-1]);

        gauss = random.boxMuller()[0]; // second gaussian law independant of the one used for V
        simulationAsset[i] = S(gauss, simulationAsset[i-1], simulationVariance[i-1], simulationVariance[i]);
    }
    
    return simulationAsset;
}

double TGScheme::V(double gauss, double v)
{
    double *tab = TG_moments(v);
    return max(tab[0] + tab[1] * gauss, 0.0);
}

double* TGScheme::TG_moments(double v)
{
    double *moment = chi_moments(v);
    double phi = moment[1] / moment[0];
    double r = newton_search(phi);
    double *res = new double [2];
    
    if(r == pow(10,10000))
    {
        res[0] = sqrt(moment[0]);
        res[0] = 0.0;
        return res;
    }
    else
    {
    res[0] = (r * sqrt(moment[0])) / (density(r) + r * CDF(r));
    res[1] = ((1 / sqrt(phi) )* sqrt(moment[1])) / (density(r) + r * CDF(r));
    return res;
    }
}

double TGScheme::newton_search(double phi)
{
    double tol = pow(10,-10);
    double x_current;
    double x_next = 10; // the zero of the equation seems to be proche to 10 according to prior simulations.
    double error;
    int it = 0;
    int max_it = 100000; //////////////ITERATION //////////////////////
    
    do
    {
        x_current = x_next;
        x_next = x_current - equation(x_current, phi) / equation_derive(x_current, phi);
        error = abs(x_next - x_current);
        it++;
    }
    while (error > tol && it < max_it && isnan(x_next) == 0);
    
    return x_current;
}

double TGScheme::equation(double r, double phi)
{
    return r * density(r) + CDF(r) * (1 + pow(r,2)) - (1 + phi) * pow( density(r) + r*CDF(r), 2 );

}

double TGScheme::equation_derive(double r, double phi) // analytic used but approximate derivative can also be used
{
    return density(r) + r * density_derive(r) + density(r) * (1 + pow(r, 2)) + CDF(r) * 2 * r - (1 + phi) * 2 *(density(r) + r*CDF(r)) * (density_derive(r) + CDF(r) + r * density(r));
    //return ( equation(r + 0.000001, phi) - equation(r, phi) ) / 0.000001;
}

double TGScheme::density(double x)
{
    return (1/sqrt(2 * M_PI)) * exp(-pow(x,2) / 2);
}

double TGScheme::density_derive(double x)
{
    return -x * (1/sqrt(2 * M_PI)) * exp(-pow(x,2) / 2);
}

double TGScheme::CDF(double x) // compute efficiently the cumulative distribution function of a normal law (see G. West: Better approximations to cumulative normal functions)
{
    static const double RT2PI = sqrt(4.0*acos(0.0));
    
    static const double SPLIT = 7.07106781186547;
    
    static const double N0 = 220.206867912376;
    static const double N1 = 221.213596169931;
    static const double N2 = 112.079291497871;
    static const double N3 = 33.912866078383;
    static const double N4 = 6.37396220353165;
    static const double N5 = 0.700383064443688;
    static const double N6 = 3.52624965998911e-02;
    static const double M0 = 440.413735824752;
    static const double M1 = 793.826512519948;
    static const double M2 = 637.333633378831;
    static const double M3 = 296.564248779674;
    static const double M4 = 86.7807322029461;
    static const double M5 = 16.064177579207;
    static const double M6 = 1.75566716318264;
    static const double M7 = 8.83883476483184e-02;
    
    const double z = fabs(x);
    double c = 0.0;
    
    if(z<=37.0)
    {
        const double e = exp(-z*z/2.0);
        if(z<SPLIT)
        {
            const double n = (((((N6*z + N5)*z + N4)*z + N3)*z + N2)*z + N1)*z + N0;
            const double d = ((((((M7*z + M6)*z + M5)*z + M4)*z + M3)*z + M2)*z + M1)*z + M0;
            c = e*n/d;
        }
        else
        {
            const double f = z + 1.0/(z + 2.0/(z + 3.0/(z + 4.0/(z + 13.0/20.0))));
            c = e/(RT2PI*f);
        }
    }
    return x<=0.0 ? c : 1-c;
    
}





QEScheme::QEScheme(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO, double threshold): HestonModel(T, N, kappa, rho, sigma, theta, r), _sO(sO),  _vO(vO), _threshold(threshold)
{
    
}

double* QEScheme::log_asset_simulation()
{
    double gauss;
    double * simulationAsset = new double [_N+1];
    double * simulationVariance = new double [_N+1];
    simulationAsset[0] = log(_sO); // BE CAREFULL : we take the log initial price
    simulationVariance[0] = _vO;
    
    for (int i=1; i<=_N; i++)
    {
        gauss = random.boxMuller()[0]; //first gaussian law
        simulationVariance[i] = V(gauss, simulationVariance[i-1]);
        
        gauss = random.boxMuller()[0]; // second gaussian law independant of the one used for V
        simulationAsset[i] = S(gauss, simulationAsset[i-1], simulationVariance[i-1], simulationVariance[i]);
    }
    
    return simulationAsset;
}

// V est calculé en fonction de la poisition de phi en fonction de phi_critics (appelé threshold)
double QEScheme::V(double gauss, double v)
{
    
    double *moment = chi_moments(v);
    double phi = moment[1] / moment[0];
    
    double U = random.randomUniform(); // uniform
    
    if (phi <= _threshold)
    {
        return next_below(sqrt(moment[0]), gauss, phi);
    }
    else
    {
        return next_above(sqrt(moment[0]), U, phi);
    }
    
}

// Fonction Next_below qui suit la procédure de calcul de Next_V lorsque phi<=phi_critic
double QEScheme::next_below(double m, double gauss, double phi)
{
    double *tab = compute_a_b(m, phi);
    double b = tab[0];
    double a = tab[1];
    
    return a * pow(b + gauss, 2);
    
}

// Fonction calculant a et b dans le cadre de la fonction Next_below
double* QEScheme::compute_a_b(double m, double phi)
{
    double *res = new double[2];
    double tempo = 2.0 / phi;
    
    res[0] = sqrt(tempo - 1.0 + sqrt(tempo) * sqrt(tempo - 1)); // b
    res[1] = m / (1.0 + pow(res[0],2)); // a
    
    return res;
}

// Fonction Next_above qui suit la procédure de calcul de Next_V lorsque phi>phi_critic
double QEScheme::next_above(double m, double U, double phi)
{
    double *tab = compute_p_beta(m, phi);
    double p = tab[0];
    double beta = tab[1];
    
    return inverse_cumulative(U, p, beta);
}

double* QEScheme::compute_p_beta(double m, double phi) {
    
    double *res = new double[2];
    
    res[0] = (phi - 1.0) / (phi + 1.0); // p
    res[1] = 2.0 / (m * (phi + 1.0)); // beta
    
    return res;
}

double QEScheme::inverse_cumulative(double U, double p, double beta) {
    
    double res(0.0); 
    if ((U >= 0) && (U <= p)) {
        res = 0.0;
    }
    else {
        res = (1.0 / beta) * log((1 - p) / (1 - U)); 
    }
    
    return res;
}
