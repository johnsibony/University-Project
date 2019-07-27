#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;

void ecrireNb(ofstream & fichier) {  // EXO 1 //
    int n;
    cout << "Entrez un nombre"<<endl;
    cin >> n;
    if(fichier) {
    for (int i=0; i<=n; i++) {
         fichier << i<<endl;
    }
     }
     fichier.close();
}

void lireNb(ifstream& fichierLire) {
    string lecture;
    cout << "Les chiffres du fichiers sont :"<<endl;
    while(getline(fichierLire, lecture)){
        cout << lecture << endl;
    }
}

int tailleFichier(ifstream& fichierLire){  //conaitre le nombre de ligne d'un fichier
    string lecture;
    fichierLire.clear();
    fichierLire.seekg(0, ios::beg);
    int i=0;
    while(getline(fichierLire, lecture)){
        i++;
    }
    return i;
}

void afficherNb (ifstream& fichierLire){
    string lecture;
    int i=0;  //compteur position curseur dans le fichier
    int n;
    cout << "Entrez un entier inferieur à " << tailleFichier(fichierLire) << endl;
    cin >> n;
    fichierLire.clear();
    fichierLire.seekg(0, ios::beg);
    while(getline(fichierLire, lecture)){
        i++;
        if(i==n){
            cout << lecture<< endl;
        }
    }
}

void assombrirEclaircir(ofstream& fichierEcrire, ifstream& fichierLire, double a){  //EXO 2 avec a=0.5 ou 2 selon assombir ou eclaircir
    string lecture;
    int longueur;
    int hauteur;
    double maximum;
    fichierLire >> lecture;
    fichierLire >> longueur;
    fichierLire >> hauteur;
    fichierLire >> maximum;
    fichierEcrire <<lecture <<endl <<longueur <<" "<<hauteur<<endl<<maximum<<endl;
    int rouge;
    int vert;
    int bleu;
    for(int i=1; i<=hauteur; i++){
        for(int j=1; j<=longueur; j++){
            fichierLire >> rouge;
            fichierLire >> vert;
            fichierLire >> bleu;
            fichierEcrire << (int)min(rouge*a,maximum)<< " "<< (int)min(vert*a,maximum)<< " "<< (int)min(bleu*a,maximum) << " "<< " "<< " ";
        }
        fichierEcrire << endl;
    }
}

void flou(ofstream& fichierEcrire, ifstream& fichierLire) { //technique : créer un cadre de pixel egale à 0 qui entoure l'image pour ne plus avoir de probleme de definition sur les bords (il y aura toujours 9voisins)
    string lecture;
    int longueur;
    int hauteur;
    double maximum;
    fichierLire >> lecture;      //on ecrit l'en tête sur le nouveau fichier
    fichierLire >> longueur;
    fichierLire >> hauteur;
    fichierLire >> maximum;
    fichierEcrire <<lecture <<endl <<longueur <<" "<<hauteur<<endl<<maximum<<endl;
    int rouge;
    double rougeMoyenne=0;       //nouvelles valeures floutées des pixels
    int vert;
    double vertMoyenne=0;
    int bleu;
    double bleuMoyenne=0;
    double tab[(hauteur+2)*(longueur+2)][3]; // +2 car on rajoute un "cadre/contour" à l'image qui vaut 0 partout
    for(int i=0; i<(hauteur+2)*(longueur+2); i++){ //on remplit le tableu RVB des pixels
        if(i%(longueur+2)==0 || (i+1)%(longueur+2)==0||i<longueur+2||i>=(longueur+2)*(hauteur+1)){ //contour vaut 0
            tab[i][0]=0;
            tab[i][1]=0;
            tab[i][2]=0;

        }
        else{
        fichierLire >> rouge;
        tab[i][0]=rouge;
        fichierLire >> vert;
        tab[i][1]=vert;
        fichierLire >> bleu;
        tab[i][2]=bleu;
            
        }
    }
    
    for(int i=0; i<(hauteur+2)*(longueur+2); i++){
        if(i%(longueur+2)!=0 && (i+1)%(longueur+2)!=0 && i>=longueur+2 && i<(longueur+2)*(hauteur+1)){ //on est pas dans le cadre
            rougeMoyenne=(tab[i-(longueur+2)-1][0]+tab[i-(longueur+2)][0]+tab[i-(longueur+2)+1][0]+tab[i-1][0]+tab[i][0]+tab[i+1][0]+tab[i+(longueur+2)-1][0]+tab[i+(longueur+2)][0]+tab[i+(longueur+2)+1][0])/9;//moyenne 9pixels voisins
            vertMoyenne=(tab[i-(longueur+2)-1][1]+tab[i-(longueur+2)][1]+tab[i-(longueur+2)+1][1]+tab[i-1][1]+tab[i][1]+tab[i+1][1]+tab[i+(longueur+2)-1][1]+tab[i+(longueur+2)][1]+tab[i+(longueur+2)+1][1])/9;
            bleuMoyenne=(tab[i-(longueur+2)-1][2]+tab[i-(longueur+2)][2]+tab[i-(longueur+2)+1][2]+tab[i-1][2]+tab[i][2]+tab[i+1][2]+tab[i+(longueur+2)-1][2]+tab[i+(longueur+2)][2]+tab[i+(longueur+2)+1][2])/9;
            fichierEcrire << (int)rougeMoyenne << " " << (int)vertMoyenne << " " << (int)bleuMoyenne <<endl;

    }
    
}
}

    int main () {
        ofstream fichierEcrire("/Users/johnsibony/desktop/C++/TP5/TP5/fichier.txt");
        ifstream fichierLire("/Users/johnsibony/desktop/C++/TP5/TP5/fichier.txt");
        ecrireNb(fichierEcrire);
        lireNb(fichierLire);
        afficherNb (fichierLire);
        ifstream fichierChat("/Users/johnsibony/desktop/C++/TP5/TP5/cat-blur.txt"); //fichier example chat converti en .txt
        ofstream fichierAssombrir("/Users/johnsibony/desktop/C++/TP5/TP5/newppmassombrir.txt");
        ofstream fichierEclaircir("/Users/johnsibony/desktop/C++/TP5/TP5/newppmeclaircir.txt");
        ofstream fichierFlou("/Users/johnsibony/desktop/C++/TP5/TP5/newppmflou.txt");
        assombrirEclaircir(fichierAssombrir,fichierChat,0.5);
        fichierChat.clear();
        fichierChat.seekg(0, ios::beg);
        assombrirEclaircir(fichierEclaircir,fichierChat,2);
        fichierChat.clear();
        fichierChat.seekg(0, ios::beg);
        flou(fichierFlou,fichierChat);
    }

