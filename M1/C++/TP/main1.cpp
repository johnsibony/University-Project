#include <iostream>
#include <iomanip>

using namespace std;

#include <stdlib.h>

#include <time.h>

#include <math.h>

int factorielle (const int n) {
    
    int a=1;
    
    for (int i=0; i<n; i++) {
        
        a=a*(i+1);
        
    }
    
    return a;
    
}


void devinette () {
    
    srand(time(NULL));
    
    int a = rand() % 21;      // générer un nb entre 0 et 20, rand compris entre 0 inclu et RAND_MAX exclu  //
    
    int n;                    // nb entré par l'utilisateur //
    
    int i=1;                  // compteur du  nombre de coup nécéssaire //
    
    cout << "Devinez le nombre compris entre 0 et 20 inclus" << endl;
    
    cin >> n;
    
    while (n !=a) {
        
        i++;
        
        cout << "Ce n'est pas le nombre cherché" << endl;
        
        if (n<a) {
            
            cout << "Le nombre est plus grand" << endl;
            
            cin >> n;
            
        }
        
        if (n>a) {
            
            cout << "Le nombre est plus petit" << endl;
            
            cin >> n;
            
        }
        
    }
    
    cout << " Bravo, vous avez gagné en " << i << " couts!" << endl;
    
}

bool nbPremier(const int n) {
    
    if (n==1) {
        
        return false;
        
    }
    if (n==2) {
        
        return true;
        
    }
    
    for (int i=2; i<=(int)sqrt(n); i++) {
        
        if (n%i==0) {
            
            return false;
            
        }
        
    }
    
    return true;
    
}


double approxPi() {
    
    double pi=0; // estimateur de Pi par Monte Carlo //
    
    const int n=100000000; // grand nombre de simulation iid //
    
    srand(time(NULL));
    
    for (int i=0; i<=n; i++) {
        
        double u= (double(rand())/RAND_MAX)*2;   // rand_max dep de XCODE mais tres grand donc
        
        double v= (double(rand())/(double)RAND_MAX)*2; // Simulation de couple de loi uniforme iid sur (0,2) //
        
        if (sqrt(pow((u-1),2) + pow((v-1),2))<=1) { // Si le couple appartient à la boule de centre 1 et de rayon 1 //
            
            pi++;                           // on compte le nombre de couple appartenant à la boule //
            
        }
        
    }
    
    return 4*pi/n;
    
}



int main(){
    
    int n;
    
    cout << "Donnez un entier naturel" << endl;
    
    cin >> n;
    
    cout << "Le factorielle de votre nombre vaut :" << factorielle (n) << endl;
    
    devinette();
    
    cout << "Donnez un entier naturel" << endl;
    
    cin >> n;
    
    cout << "Votre nombre est premier? 1 si oui, 0 si non" << endl << nbPremier(n) << endl;
    
    cout << "l'approximation de Pi par la méthode de Monte Carlo vaut : " << setprecision(9) << approxPi() << endl;
    
}
