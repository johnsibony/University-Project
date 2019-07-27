
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tâche 7 à 9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nombre des noeuds horizontaux de la grille
N=16;

% NOmbre des noeuds verticaux de la grille
M=8;

% Durée 
T=10;


% Pas spatial
h=4/N;

% Valeurs initiales

% Température initiale
u0=290;
nu=10;
g1=28.5;
g2=-28.5;
g3=0;



% Les quatres listes de numérotations Nord/ Sud/Est et Ouest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%liste Nord
Nord=[N+2:(N+1)*(M+1) zeros(1,N+1)];

%Liste Sud
Sud=[zeros(1,N+1) 1:N*M+M];

% LIste Ouest 
Ouest=[0 1:(N+1)*(M+1)-1];
for k=1:M
j=(N+1)*k+1;
Ouest(j)=0;
endfor

%Liste Est
Est=[2:(N+1)*(M+1) 0];
for k=1:M+1
j=(N+1)*k;
Est(j)=0;
endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction de La matrice A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Le domaine Omega
%%%%%%%%%%%%%%%%%%
A=diag(4/h^2*ones(1,(N+1)*(M+1)));
A=sparse(A);
for k=1:(N+1)*(M+1)
if Nord(k)~=0
A(k,Nord(k))=-1/h^2;
endif

if Sud(k)~=0
A(k,Sud(k))=-1/h^2;
endif

if Est(k)~=0
A(k,Est(k))=-1/h^2;
endif

if Ouest(k)~=0
A(k,Ouest(k))=-1/h^2;
endif

endfor

%% Les 4 coin de la grille
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A(1,Est(1))=-2/h^2;
A(1,Nord(1))=-2/h^2;


A(N+1,Ouest(N+1))=-2/h^2;
A(N+1,Nord(N+1))=-2/h^2;


A(N*M+M+1,Sud(N*M+M+1))=-2/h^2;
A(N*M+M+1,Est(N*M+M+1))=-2/h^2;


A((N+1)*(M+1),Sud((N+1)*(M+1)))=-2/h^2;
A((N+1)*(M+1),Ouest((N+1)*(M+1)))=-2/h^2;


%% Les points appartenant à Gamma 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:M+1
j=(N+1)*k;
A(j,Ouest(j))=-2/h^2;
endfor


for k=1:M
j=(N+1)*k+1;
A(j,Est(j))=-2/h^2;
endfor

%Les points appartenant à Gamma 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=N*M+M+2:(N+1)*(M+1)-1
A(k,Sud(k))=-2/h^2;
endfor


%Les points appartenant à Gamma 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=2:N
A(k,Nord(k))=-2/h^2;
endfor


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Construction du vecteur F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialisation de F
F=zeros((N+1)*(M+1),1);

% Les quatres coins
F(1)=-2*nu*(g2+g3)/h;
F(N+1)=2*nu*(-g2+g3)/h;
F(N*M+M+1)=2*nu*(g1-g3)/h;
F((N+1)*(M+1))=2*nu*(g1+g3)/h;

%Les points appartenant à Gamma 1

for k=N*M+M+2:(N+1)*(M+1)-1
F(k)=2*nu*(g1)/h;
endfor

%Les points appartenant à Gamma 2

for k=2:N
F(k)=-2*nu*(g2)/h;
endfor

% Les points appartenant à Gamma 3

% à droite de la grille
for k=2:M
j=(N+1)*k;
F(j)=2*nu*(g3)/h;
endfor

% à gauche de la grille
for k=1:M-1
j=(N+1)*k+1;
F(j)=-2*nu*(g3)/h;
endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Implémentation du  theta shéma 

%Initialisation de u
u=ones((N+1)*(M+1),1)*u0;

%Valeur de theta
theta=1;

% Pas temporel
dt=0.1;

% Nombre d'itéraiton temporel
kk=floor(T/dt)-1;

% Vérification de la stabilité
if theta<1/2
[V,D]=eigs(A);
Vdisc=sort(diag(D)');
lambda_max=max(Vdisc);
stability=2/(nu*lambda_max*(1-2*theta));
if dt>stability
error('methode non stable, diminuez dt')
endif
endif

%theta shéma
w=u;
w_mean=mean(u);
C=inv(eye((N+1)*(M+1))+dt*nu*theta*A)*(eye((N+1)*(M+1))-dt*nu*(1-theta)*A);
D=inv(eye((N+1)*(M+1))+dt*nu*theta*A);
for k=1:kk
u=C*u+dt*D*F;
w=[w,u];
w_mean=[w_mean,mean(u)];
endfor


















