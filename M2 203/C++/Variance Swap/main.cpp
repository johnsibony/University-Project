#include <iostream>
#include <time.h>
#include "varianceSwapPricer.hpp"
#include "europeanOptionPricer.hpp"

using namespace std;



int main()

{
    srand((unsigned int)time(NULL)); // initialisation of the rand generator
    string res;
    
    cout << "Would you like to price a Variance swap or European option ? (swap, option)" << endl;
    cin >> res;
    
    if(res == "swap")
    {
    
    VarianceSwapPricer swap;
        
        /* You can use 2 others constructors :
         
         * VarianceSwapPricer(double T, int N, double sO, double vO)
         * VarianceSwapPricer(double T, int N, double kappa, double sigma, double theta, double r, double rho, double sO, double vO, double threshold)
         
         with sO, vO and threshold the initial asset price, initial variance price and the critical value of Phi in QE scheme
         
         */
    try
    {
        complex<double> price = swap.get_analytic_price("re"); // choose "re" or "im" or no argument for respectively real or imaginary or (total) complexe part.
        cout << "Analytic Strike Price = " << price << endl ;
    }
    catch(string const& error)
    {
        cerr << error << "in Analytic method" << endl;
    }
        cout << "What kind of scheme would you use ? (euler, tg, qe)" << endl;
        cin >> res;
        
        try
        {
            double price = swap.get_montecarlo_price(res, 500000); // number of MC simulation
            cout << res << " Scheme Monte Carlo Strike Price = " << price << endl ;
        }
        catch(string const& error)
        {
            cerr << error << " in '" << res << "' scheme" << endl;
        }
    }

    else if(res == "option")
    {
        
        EuropeanOptionPricer option;
        
        /* You can use 2 others constructors :
         
         * EuropeanOptionPricer(double T, int nbMC, double sO, double vO, double K)
         * EuropeanOptionPricer(double T, int nbMC, double kappa, double sigma, double theta, double r, double rho, double sO, double vO, double K, double threshold)
         
         with sO, vO, K, threshold and nbMC the initial asset price, initial variance price, strike price, critical value of Phi in QE scheme and number of MonteCarlo simulation
         
         */
        try
        {
            complex<double> price = option.get_analytic_price(1000000); // number of rimann integral simulation
            cout << "Analytic Strike Price = " << price << endl ;
        }
        catch(string const& error)
        {
            cerr << error << "in Analytic method" << endl;
        }
        
        cout << "What kind of scheme would you use ? (euler, tg, qe)" << endl;
        cin >> res;
        
        try
        {
            double price = option.get_montecarlo_price(res, 500000); // number of MC simulation
            cout << res << " Scheme Monte Carlo Strike Price = " << price << endl ;
        }
        catch(string const& error)
        {
            cerr << error << " in '" << res << "' scheme" << endl;
        }
    }
    else
    {
        cout << "Invalid answer" << endl;
    }
    

}
