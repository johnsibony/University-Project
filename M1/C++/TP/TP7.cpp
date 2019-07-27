/// TD7

#include <iostream>
#include <cstdlib>

using namespace std;


struct Node{
    int data;
    struct Node* next,*prev;
};


class BigNum{
    struct Node *head,*tail; // liste chainÈe bidirectionnellle
    public :
    BigNum(int); // constructeur avec int en param
    BigNum(BigNum&); // constructeur avec pointeur sur BigNum en param
    BigNum& operator=(BigNum); //overloading de loperateur egale
    ~BigNum(); // destructeur
    friend BigNum& operator+(BigNum, BigNum);
    friend BigNum& operator*(BigNum, BigNum);
    friend BigNum& operator*(BigNum, int);
    friend ostream& operator <<(ostream &, BigNum); //overloading affichage d'un BigNum
    
    
};


BigNum::BigNum(int n){
    int a = n;
    int r = a%10;
    this->tail = new Node;
    this->head = this->tail;
    this->tail->next = this->tail->prev = NULL;
    Node* p = this->head;
    Node*tmp;
    p->data = r;
    while((a=a/10)!=0){
        r = a%10;
        tmp =  p;
        p->prev = new Node;
        p->prev->data = r;
        p = p->prev;
        p->next = tmp;
        p->prev = NULL;
    }
    
    this->head = p;
    
}


BigNum::BigNum(BigNum& B){ /// copy constructor
    Node* b = B.head ;
    this->head = new Node;
    this->head->prev = this->head->next = NULL;
    Node* a = this->head;
    Node* tmp;
    a->data = b->data;
    b = b->next;
    
    while(b!=NULL){
        a->next = new Node;
        tmp = a;
        a = a->next;
        a->prev = tmp;
        a->data = b->data;
        
        b = b->next;
    }
    
    this->tail = a;
    this->tail->next = NULL;
    
}

BigNum& BigNum::operator=(BigNum A){
    //        cout<<"on est dans loverloading ="<<endl;
    Node* a = A.tail;
    this->tail = new Node;
    this->head = this->tail;
    this->head->prev = this->head->next = NULL;
    Node*p = this->head;
    Node*tmp;
    p->data = a->data;
    //        cout<<"hello 2"<<endl;
    a = a->prev;
    //        cout<<"hello1"<<endl;
    while(a!=NULL){
        p->prev = new Node;
        tmp = p;
        p=p->prev;
        p->data = a->data;
        p->next = tmp;
        p->prev = NULL;
        a = a->prev;
    }
    this->head = p;
    return *this;
}

BigNum& operator+(BigNum A, BigNum B){
    BigNum* C = new BigNum(0);
    Node * c = C->tail;
    
    Node * a = A.tail;
    Node * b = B.tail;
    Node*tmp;
    int s,r = 0;
    /// r est indicatrice de retenue
    
    s = a->data + b->data; // on suppose que les deux bignum en parametre comporte au moins un chiffre
    if(s>9){
        r++;
        s=s%10;
    }
    c->data = s;
    a = a->prev;
    b = b->prev;
    
    while(a!= NULL || b!=NULL){ // la boucle sarrete quand on a parcouru les deux listes en entier
        tmp = c;
        c->prev = new Node;
        c = c->prev;
        c->next = tmp;
        c->prev = NULL;
        
        if(a!=NULL){
            if(b!=NULL){s = a->data + b->data;}
            else{s=a->data;}
        }
        else if(a==NULL){
            s = b->data;
        }
        
        
        if(r != 0){ // il y a une retenue
            s = s+r;
            r = 0;
            
        }
        if(s>9){
            r = s/10;
            s=s%10;
        }
        
        c->data = s;
        
        if(a!=NULL){a = a->prev;}
        if(b!=NULL){b = b->prev;}
    }
    
    if(r!=0){
        tmp = c ;
        c->prev = new Node;
        c = c->prev;
        c->data = r;
        c->next = tmp;
        c->prev = NULL;
    }
    C->head = c;
    return *C;
}

BigNum& operator*(BigNum A,int n){
    //    cout<<"on est dans loverloading de *"<<endl;
    //cout<<"     a->data ="<<n<<endl;
    BigNum* B = new BigNum(0);
    for(int i=1;i<=n;i++){
        *B = (A+*B);
        
    }
    //    cout<<"        *B:"<<*B<<endl;
    return *B;
}

BigNum& operator*(BigNum a,BigNum b){
    BigNum* result = new BigNum(0);
    BigNum* subresult=new BigNum(0);
    Node*tmp;
    tmp=a.tail;
    int i=0;
    while(tmp){
        *subresult=b*(tmp->data);
        for (int j=0;j<i;j++){
            *subresult=*subresult*10;
            
        }
        *result=*result+*subresult;
        i=i+1;
        tmp=tmp->prev;
        
    }
    return *result;
    
    
}

ostream& operator<<(ostream& sortie,BigNum A){
    //    out<<"on est dans l'overloading de << "<<endl;
    Node* a = A.head;
    while(a!=NULL){
        sortie << (a->data);
        a = a->next;
    }
    return sortie;
}


BigNum::~BigNum(){
    Node* p = this->head;
    Node*prec;
    while(p!=NULL){
        prec = p;
        p = p->next;
        delete(prec);
    }
}




BigNum& fact(int n){
    BigNum* result = new BigNum(1);
    if(n==0 || n==1){
        return *result;
    }
    for(int i=2;i<=n;i++){
        *result = *result*i;
    }
    return *result;
}


BigNum& fibo(int n){
    BigNum *a = new BigNum(0);
    BigNum *b = new BigNum(1);
    BigNum *temp = new BigNum(0);
    for(int i=2;i<=n;i++){
        *temp = *b;
        //        cout<<"t"<<*temp<<endl;
        *b = *b+*a;
        *a = *temp;
    }
    return *b;
}


int main(){
    //    int x =4;
    //    int *p = &x;
    //    cout<<p<<endl;
    BigNum F(123456789);
    BigNum G(987654321);
    
    cout<<"F = "<<F<<" G = "<<G<<endl;
    
    cout<<"la somme de F et G vaut : "<<F+G<<endl;
    //        cout<<"la multiplication de F et 5 vaut : "<<F*78<<endl;
    
    cout<<"la multiplication de F et G vaut : "<<F*G<<endl;
    
    
    
    BigNum B(fact(30));
    
    cout<<"facto de 30 au carre vaut : "<<B*B<<endl;
    
    BigNum C(fibo(60));
    cout<<"fibo de 60 vaut:"<<C<<endl;
    
    //        cout<<"factoriel de trente"<<result<<endl;
    return 0;
}
