/*
 Ce fichier contient 2 classes utiles à la poursuite du projet. L'une permettant d'élargir les propriétés de la class complexe, l'autre de créer de l'aléatoire
 */

#ifndef utilsClass_hpp
#define utilsClass_hpp

#include <complex>




/* 
 New complex class allowing classical operation between int/double and complex object.
 The user will only have to create a complexe number i (0 + 1*i) and the class will 
 automatically recognize its real part and imaginary part.
 
 Example : Instead of creating the complexe 1+2i, we can create the number i of the class My_complexe
 and then transform it through 1+2i
*/

class My_complex: public std::complex<double>
{

    public:
    
    My_complex(double a, double b);
    My_complex(complex<double> a);

    private:
    
    /*
     Usual operand between double and complex
     The friend concept allow us to use private atribute in the function outside of the class
    */
    friend My_complex operator+ (const double, const My_complex);
    friend My_complex operator- (const double, const My_complex);
    friend My_complex operator* (const double, const My_complex);
    
    /*
     Usual operand between int and complex
     The friend concept allow us to use private atribute in the function outside of the class
    */
    friend My_complex operator+ (const int, const My_complex);
    friend My_complex operator- (const int, const My_complex);
    friend My_complex operator* (const int, const My_complex);

};

/*
 We define the function charged to the complex operations
*/
My_complex operator* (const double nb, const My_complex z);
My_complex operator+ (const double nb, const My_complex z);
My_complex operator- (const double nb, const My_complex z);

My_complex operator* (const int nb, const My_complex z);
My_complex operator+ (const int nb, const My_complex z);
My_complex operator- (const int nb, const My_complex z);


/*
 A random class to simulate gaussian couple with correlation with the Box Muller method
*/
class RandomTool
{
    
    public:
    
    RandomTool();
    RandomTool(double rho);
    double randomUniform(); // simulation of the uniform law
    double* boxMuller();
    
    private:
    
    double _rho; // correlation coeficient
    
};

#endif
