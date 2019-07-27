#include <iostream>
using namespace std;


class BigNum {
    
private:
    int *tab;   //l'unique attribut est le tableau ou est stcoké le bigNum
    
public:
    BigNum() {   //constructeur par defaut
        tab=new int[100];
        for (int i=99; i>=0; i--){
            tab[i]=0;
        }
    }
    BigNum(int a) {   //constructeur avec argument
        tab=new int[100];
        for (int i=99; i>=0; i--){
            tab[i]=a%10;
            a=a/10;
        }
    }
    friend BigNum operator+(BigNum, BigNum );
    friend BigNum operator*(BigNum const&, BigNum const&);
    friend BigNum operator*(BigNum const&, int const&);
    friend bool operator<(BigNum, BigNum);
    friend ostream& operator<<(ostream&, BigNum const&);
};

BigNum operator+(BigNum a, BigNum b){
    BigNum resultat;
    int unite=0;
    for (int i=99; i>=0; i--){
        if((a.tab[i]+b.tab[i]+unite)>9){  //si la somme est composé de 2chiffres
            resultat.tab[i]=(a.tab[i]+b.tab[i]+unite)%10;  //on prend l'unité du nombre à 2chiffres
            unite=1; //on rajoute une unité pour la prochiane somme
        }
        else{
            resultat.tab[i]=a.tab[i]+b.tab[i]+unite;  //sinon on ecrit l'unique chiffre
            unite=0; //et on reinitialise l'unité à 0
        }
    }
    return resultat;
}


BigNum operator*(BigNum const&a, BigNum const&b){  //rem : on doit avoir (nombre de chiffre de a + nombre de chiffre de b)<100
    BigNum resultat;
    for (int i=99; i>=0; i--){  //les chiffres de b
        int unite=0; // la retenu de la multiplication
        BigNum resultatintermediaire;
        for (int j=99; j>=99-i; j--){   //les chiffres de a
            if((a.tab[j]*b.tab[i])+unite>9){ //si le produit +la retenue est composée de 2chiffres
                resultatintermediaire.tab[j-(99-i)]=((a.tab[j]*b.tab[i])+unite)%10; //on ecrit l'unité du nombre à 2chiffres. On remplit le tableau de la case 0 à 99-i car on laisse les 0 de 99-i+1 à 99
                unite=((a.tab[j]*b.tab[i])+unite)/10;  //la retenue est la dizaine
            }
            else {
                resultatintermediaire.tab[j-(99-i)]=(a.tab[j]*b.tab[i])+unite; //on ecrit ce chiffre unique avec la retenue
                unite=0;  //sinon l'unité vaudra 0 car composé de 1chiffre unique
            }
        }
        for (int j=99; j>=0; j--){ //puis il faut additionner chaque multiplication entres elles
            if((resultat.tab[j]+resultatintermediaire.tab[j])>9){
                (resultatintermediaire.tab[j-1])+=1;
                resultat.tab[j]=(resultat.tab[j]+resultatintermediaire.tab[j])%10;
            }
            else {
                resultat.tab[j]+=resultatintermediaire.tab[j];
            }
        }
    }
    return resultat;
}

BigNum operator*(BigNum const&a, int const&b ){
    BigNum resultat(b);
    resultat=resultat*a;
    return resultat;
}

bool operator<(BigNum a, BigNum b){
    int carda=0;  //nombre de 0 avant le chiffre de a
    int cardb=0;//nombre de 0 avant le chiffre de b
    while(a.tab[carda]==0){
        carda++;
    }
    
    while(b.tab[cardb]==0){  //refait pareil pour b
        cardb++;
    }
    if (carda>cardb){
        return true; //le tableau a à plus de 0 avant donc moins de chiffre donc plis petit que b forcement
    }
    if(carda<cardb) { //et inversement
        return false;
    }
    for (int i=carda; i<=99; i++){ //on commence a partir du premier chiffre en commun
        if (a.tab[i]<b.tab[i]){ //si le premier chiffre (en partant a gauche) de a est plus petit c'est bon
            return true;
        }
        if (a.tab[i]>b.tab[i]){  //sinon c'est faux
            return false;
        }
        // si c'est egale on regarde le chiffre suivant avec la boucle for
    }
    return false; //si il sont exactement egaux alrs c'est faux
}

ostream& operator<<(ostream& fichier, BigNum const& a){
    int carda=0;  //nombre de 0 avant le chiffre de a
    while(a.tab[carda]==0){
        carda++;
    }
    for (int i=carda; i<=99; i++){
        fichier <<a.tab[i];
    }
    return fichier;
}

BigNum factorielle(int a){
    BigNum big(1);
    if(a==0){
        //cout <<big<<endl;
        return big;
    }
    
    for (int i=1; i<=a;i++){
        big=big*i;
        //cout <<big<<endl;
    }
    return big;
}

BigNum fibonacci(int a){
    BigNum big1(0);
    BigNum big2(1);
    BigNum bigtotal;
    if(a==1){
        return BigNum(0);
    }
    if(a==2){
        return BigNum(1);
    }
    for (int i=1; i<=a-2;i++){
        bigtotal=big1+big2;
        BigNum temp=big2;
        big2=bigtotal;
        big1=temp;
    }
    return bigtotal;
}

int main() {  //On donne des exemples pour aider la prof
    cout << "Entrez 2 entiers" << endl;
    int a,b;
    cin >> a>>b;
    BigNum numero1(a);
    BigNum numero2(b);
    cout << "Premier BigNum: "<<numero1<<endl<< "Second BigNum: "<<numero2<<endl<< "Somme des 2 BigNum: "<<numero1+numero2<<endl<< "Multiplication des 2 BigNum: "<<numero1*numero2<<endl<< "Multiplication BigNum et entier: "<<numero1*b<<endl<< "Premier Bignum inferieur au second? 1=oui/0=non: "<<(numero1<numero2)<<endl<<endl<<endl<<"Le factoriel de 30 vaut : "<<factorielle(30)<<endl<<endl<<endl<<"Le 60ème nombre de fibonacci (sachant que le premier terme et le seconde valent 0 et 1) est : "<<fibonacci(60)<<endl;
    
    return 0;
}
