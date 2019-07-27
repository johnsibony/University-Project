/*
 Ce fichier permet de simuler la variance des prix de l'asset S selon la méthode analytique
 */


#ifndef analyticPrice_hpp
#define analyticPrice_hpp

#include "utilsClass.hpp"


/*
 Class to define the function a,b and g of the analytic method
*/

class All_functions
{
    
    public:
    
    All_functions(double kappa, double rho, double sigma);
    My_complex a(double omega);
    My_complex b(double omega);
    My_complex g(double omega);
    
    protected:
    
    double _kappa;
    double _rho;
    double _sigma;
    
};

/*
 Class to define the function C with its first and second derivative
*/

class C: public All_functions
{
    
    public:  // pour que le prof puisse afficher les résultats de fonction C
    
    C(double kappa, double rho, double sigma, double theta, double r);
    My_complex value_C(double tau, double omega);
    My_complex first_derivative(double tau, double h);
    My_complex second_derivative(double tau, double h);
    
    protected:
    
    double _r;
    double _theta;
    
};

/*
 Class to define the function D with its first and second derivative
*/


class D: public All_functions
{
    
    public:  // pour que le prof puisse afficher les résultats de la fonction D
    
    D(double kappa, double rho, double sigma);
    My_complex value_D(double tau, double omega);
    My_complex first_derivative(double tau, double h);
    My_complex second_derivative(double tau, double h);
    
};

/*
 Class to define the variance of the asset on each time path
*/

class Analytic
{
    
    public:
    
    Analytic(double T, int N, double vO, double kappa, double rho, double sigma, double theta, double r);
    double time(int i); // return the ith point of the constant segmentation time
    double time_path();
    double c(int i);
    double W(int i);
    double q();
    My_complex variance_asset_initial(); // compute the initial variance asset (i=0)
    My_complex variance_asset(int i); // compute the others variance asset (i>0)
    
    private:
    
    double _T;
    int _N;
    double _vO;
    double _kappa;
    double _rho;
    double _sigma;
    double _r;
    double _theta;
    C C;
    D D;
    
};

#endif
