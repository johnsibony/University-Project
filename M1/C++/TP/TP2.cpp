#include <iostream>
#include <iomanip>
using namespace std;
#include <stdlib.h>
#include <time.h>
#include <math.h>

bool nbPremier(int const& n) {
    
    if (n==1) {
        return false;
              }
    if (n==2) {
        return true;
              }
    for (int i=1; i<=(int)sqrt(n); i++) {
        
        if (n%(i+1)==0) {
            return false;
                        }
                                        }
        return true;
                             }

bool goldbachUnite (int const& n) {   // vérifie si 1 seul entier n est de goldbach //
    if (n<=2 || n%2!=0) {
        return false;
    }
    else {
    for (int j=4; j<=n ; j+=2) {
        for (int i=2; i<=j-2; i++) {     // on s'arrête à j-2 car j et 1 ne sont pas premiers //
            if (nbPremier(i) && nbPremier(j-i)) {   // les 2 entiers doivent etre premier //
                cout << j << "=" << i << "+" << j-i << endl;
                break;
                                                } // on sort de la première boucle dès qu'on trouve un couple bon //
            if (i==j-2) return false;         // si on a parcouru toute la boucle sans break alors il n'y a pas de couple //
                                   }
                                }
    return true;                     // si on sort de la boucle prncipale alors il existe toujours un couple //
        }
                                   }
    



void renverse (int tab[],int const& size) {
    
    for (int i=0; i<size/2; i++) {      // on parcourt la première moitié du tableau //
        int a=tab[i];                   // on garde en mémoire la valeur //
        tab[i]=tab[size-i-1];           // la valeur parcourut est switcher avec sa case symétrique //
        tab[size-i-1]=a;
                                 }
    cout << "Le tableau renversé vaut :" << endl; // on affiche le tableau renversé pour aider la prof lors de sa
                                                                                                 // correction //
    for (int i=0; i<size; i++) {
        cout << tab[i];
                               }
                                          }
void itoa(int const& n, char *s){
    int a=n;                      // on copie n car n constant et on va devoir le modifier //
    char b[10]={'0','1','2','3','4','5','6','7','8','9'};
    for (int i=0; i<=4; i++) {    // 5 itération pour les 5 cases du tableau //
        *(s + 4-i)=b[a%10];          // La case i en partant de la fin prend le dernier chiffre de a //
        a=(a-a%10)/10;       // à chaque itération on enlève le dernier chiffre de a //
                             }
    for (int i=0; i<=4; i++) {   // on affiche les valeurs du tableau pour aider la prof lors de sa correction //
        cout << *(s + i);
                             }
                                }

void toBin(int const& n,int *b) {
    int a=n;                    // on copie n car n constant et on va devoir le modifier //
    for (int i=0; i<=9; i++) {   // il y a 10 chiffres pour la représentation binaire //
        *(b +9-i)=a%2;          // le reste de n sur 2 donne vaut toujours 0 si paire ou 1 si impaire //
        a=(int)a/2;             // on modife a //
                             }
    for (int i=0; i<=9; i++) {   // on affiche les valeurs du tableau pour aider la prof lors de sa correction //
        cout << *(b+i);
                             }
                          }


int main () {
    int n;
    cout << "Entrez un nombre pair (>2) pour verifier la conjecture de Goldbach. 1 si vérifié, 0 sinon " << endl;
    cin >> n;
    cout << goldbachUnite(n) << endl;
    cout << "Donnez la taille du tableau" << endl;
    cin >> n;
    int tab [n];
    cout << "Entrez les " << n << " valeurs du tableau"<< endl;
    for (int i=0; i<n; i++) {
        cin >> tab[i];
                            }
    renverse (tab, n);
    cout << endl << "Entrez un entier naturel pour afficher ses 5 derniers chiffres" << endl;
    cin >> n;
    char s;
    itoa(n, &s);
    cout << endl <<"Entrez un entier compris entre 0 et 1023 pour avoir sa representation binaire" << endl;
    cin >> n;
    int b;
    toBin(n,&b);
    cout << endl;
            }

