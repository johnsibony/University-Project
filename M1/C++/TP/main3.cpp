#include <iostream>
#include <iomanip>
using namespace std;


// Exo 1 //

int *allocationTableau(int taille) {
    int *a=new int[taille];    //allocation dynamqiue de pointeur vers tableau de taille "taille" //
    return a;                  // on fera attention à delete le pointeur plus tard //
                                   }

void remplissageTableau(int *a, int taille) {   // le tableau sera modifié //
    for (int i=0; i<taille; i++) {
        cout << "Entrez la valeur n°" << i << endl;
        cin >> a[i];
                                 }
                                        }

int sommeTableau(int *a, int taille) {
    int compteur=0;
    for (int i=0; i<taille; i++) {
        compteur += a[i];
                                 }
    return compteur;
                                        }

int *somme2Tableaux(int *u,int *v , int taille) {
    int *w=allocationTableau(taille);
    for (int i=0; i<taille; i++) {
        w[i]=u[i]+v[i];
                                 }
    return w;     // on fera attention à delete le pointeur plus tard //
                                                  }

void desallocationTableau(int *a) {
    delete [] a;
                                  }

// Exo 2 //


int **allocationTableaudb(int taille) {
    int **a=new int*[taille];
    for (int i=0; i<taille; i++) {
        a[i]=allocationTableau(taille);
                                 }
        return a;                  // on fera attention à delete le pointeur plus tard //
                                      }

void remplissageTableaudb(int **a, int taille) {
    for (int i=0; i<taille; i++) {
        cout << "Dans la ligne n°" << i << ":"<< endl;
        remplissageTableau(a[i], taille);
                                 }
                                               }

int sommeTableaudb(int **a, int taille) {
    int compteur=0;
    for (int i=0; i<taille; i++) {
        compteur += sommeTableau(a[i], taille);
                                 }
    return compteur;
                                         }

int **somme2Tableauxdb(int **u,int **v, int taille) {
    int **w=allocationTableaudb(taille);
    for (int i=0; i<taille; i++) {
        w[i]=somme2Tableaux(u[i], v[i], taille);
                                 }
    return w;     // on fera attention à delete le pointeur plus tard //
                                                     }

void desallocationTableaudb(int **a, int taille) {
    for (int i=0; i<taille; i++) {
        desallocationTableau(a[i]);
                                 }
    delete [] a;
                                                 }

// Exo 3 //

int **produitMatricielle(int **u,int **v, int taille ) { //on fait l'hypothèse que ce sont des matrice carrées de même //
                                                                // taille puisque ce n'est pas explicité dans l'énoncé //
    int **w=allocationTableaudb(taille);
    for (int i=0; i<taille; i++) {   // ligne i de la matrice W //
        for (int j=0; j<taille; j++) {  // colonne j de la matrice W //
            w[i][j]=0;
            for (int k=0; k<taille; k++) {  // indice k de la somme de la formule du produit //
               w[i][j]=w[i][j] + (u[i][k]*v[k][j]);  // formule du produit matricielle classique //
                                         }
                                     }
                                  }
    return w;
                                                        }

int main () {                 // exemple d'application de l'exo 2 et 3 avec taille=3 pour des tableaux à 2 dimensions //
    int **a=allocationTableaudb(3);
    int **b=allocationTableaudb(3);
    remplissageTableaudb(a, 3);
    remplissageTableaudb(b, 3);
    int c=sommeTableaudb(a, 3);
    cout << "La somme des éléments du premier tableau vaut : " << c << endl;
    int **d=somme2Tableauxdb(a,b,3);
    cout << "La somme des 2 tableaux vaut : " << endl;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            cout << d[i][j] << " ";
                                }
        cout << endl;
                            }
    cout << "Le produit matricielle des 2 tableaux vaut : " << endl;
    int **e=produitMatricielle(a,b,3);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            cout << e[i][j] << " ";
                                }
        cout << endl;
                            }
    desallocationTableaudb(a,3);
    desallocationTableaudb(b,3);
    desallocationTableaudb(d,3);
    desallocationTableaudb(e,3);
            }

