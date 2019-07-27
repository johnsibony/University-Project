#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

bool operateur( char *s, int j) {   //savoir si la case n°j est une opétaion (true) ou un chiffre (false) //
    if (s[j]!='0' && s[j]!='1' && s[j]!='2' && s[j]!='3' && s[j]!='4' && s[j]!='5' && s[j]!='6' && s[j]!='7' && s[j]!='8' && s[j]!='9') {
        return true;
           }
    return false;
                                 }

// fonctionnement algo exo 1 : on a une chaine de caractere au depart qu'on va séparer en 2 blocs à chaque fois par recursion. Pour cela, dès qu'on a une opération on retourne 2 blocs du tableau : 1 petit bloc à un caractère qui correspon à une des 2 opérandes et un gros bloc avec le reste qui correspond à l'autre opérande (non formé).
// Comment fabriquer ces 2 blocs?
// Si apres l'opérateur on a un chiffre alors le premier bloc sera cette opérande (petit bloc) et le second bloc correspondra à tout ce qui y'a apres ce chiffre (gros bloc)
// Si après l'opérateur on a encore un opérateur, alors le gros bloc correspondant à la premiere opérande sera fabriqué du second opérateur en question jusqu'à l'avant derniere case, et le petit bloc sera l'unique derniere case correspondant à la deuxième opérande
// Tout le problème est donc de bien savoir séparer la première opérande de la seconde
// Le bloc unique est toujours constitué d'une seule opérande et ce sera celui la qui sera affiché. Donc l'algorithme par recursion retournera toujours un tableau de taille 1 pour chaque opérande. Les opérateurs sont bien évidemment affichés juste apres que la premiere opérande est affiché.

void afficheInfix(ostream & sortie, char *s, int taille) {  //exo 1//
    if (taille==1) {
        sortie<< s[0]; //le bloc unique correspond au 1er chiffre(opérande) à afficher //
        return;   // on sort de la boucle récursive pour revenir à la précédente et afficher l'opérateur //
    }
    if (operateur(s, 1)) {   // si la case suivante est un opérateur alors la premiere opérande est le gros bloc
        if(operateur(s, taille-1)==false &&operateur(s, taille-2)==false&&operateur(s, taille-3)==true) { //on fait attention au cas ou les 3 dernière cases sont de la forme opérateur opérande1 operande 2 dans ce cas ces trois cases formes une nouvelle operande
        sortie<<"(";
        afficheInfix(sortie, s+1, taille-4);   // dans ce cas le gros bloc est le meme mais en enlevant ces trois dernieres cases
        sortie<< s[0];
        afficheInfix(sortie, s+taille-3, 3);   // et donc le petit bloc n'est plus de taille 1 mais 3
        sortie<<")";
        }
        else {  //si ce n'est pas le cas, le gros bloc est celui habitule et le peit bloc est unique (une case : la derniere ici)
            sortie<<"(";
            afficheInfix(sortie, s+1, taille-2);
            sortie<< s[0];
            afficheInfix(sortie, s+taille-1, 1);
            sortie<<")";
        }
    }
    else {  // si la case suivante est une opérande alors il correspond au petit bloc unique
        sortie<<"(";
        afficheInfix(sortie, s+1, 1);
        sortie<<s[0];
        afficheInfix(sortie, s+2, taille-2);  // le deuxième bloc (gros bloc) est le reste
        sortie<<")";
    }
    
    
}

int taille (char * s) {  //exo 3//
    int i=0;
    while (s[i]!='\0'){
        i++;
    }
    return i;
}

    int main (){
    int n;
    cout << "Entrez la taille du tableau à caractère";
    cin >>n;
    char *s=new char[n];
    cout << "Entrez les valeurs du tableau préfixe";
    for (int i=0; i<n; i++)  {
        cout << "Entrez la valeur n°" << i << endl;
        cin >> s[i];
                                  }
        cout <<"Taille du tableau : " << taille(s) << endl;;
        cout << "Tableau infixe : ";
        afficheInfix(cout,s,n);
        cout << endl;
        
    }

