clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tche 1  3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Nombre de points dans la grille
N=10;
% Pas de la grille
h=2/N;


% Les quatres listes de numrotations Nord/ Sud/Est et Ouest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Liste Nord
Nord=[N:(N-1)*(N-1) zeros(1,N-1)];

%Liste Sud
Sud=[zeros(1,N-1) 1:(N-1)*(N-2)];

%Liste Ouest
Ouest=[0 1:(N-1)*(N-1)-1];
for k=1:N-2
j=(N-1)*k+1;
Ouest(j)=0;
endfor

%Liste Est
Est=[2:(N-1)*(N-1) 0];
for k=1:N-2
j=(N-1)*k;
Est(j)=0;
endfor


% Construction de la matrice A 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Remplissage de la diagonale
A=diag(4/h^2*ones(1,(N-1)*(N-1)));

% Format creuse
A=sparse(A);

% Utilisation de la liste Nord 
for k=1:(N-1)*(N-1)
if Nord(k)~=0
A(k,Nord(k))=-1/h^2;
endif

%Utilisation de la liste Sud
if Sud(k)~=0
A(k,Sud(k))=-1/h^2;
endif

%Utilisation de la liste Est
if Est(k)~=0
A(k,Est(k))=-1/h^2;
endif

%Utilisation de la liste Ouest
if Ouest(k)~=0
A(k,Ouest(k))=-1/h^2;
endif

endfor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Une deuxime faon de construire la matrice A en utilisant l'oprateur Kron
%N=4;
%I = speye(n-1);
%e = ones(n-1,1);
%A = spdiags([e,-4*e,e],[-1,0,1],n-1,n-1);
%J = spdiags([e,e],[-1,1],n-1,n-1);
%Lh = kron(I,A) + kron(J,I);

%Etude des valeurs propres 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Les 10 plus petites valeurs propres de A
[V,D]=eigs(A,10,'sm');
% On met les valeurs prores en ordre croissant
Vdisc=sort(diag(D)');

%les valeurs propores analytiques
Vth=[pi^2/2 5/4*pi^2 5/4*pi^2 2*pi^2 13/4*pi^2 13/4*pi^2 9/2*pi^2 25/4*pi^2 25/4*pi^2 8*pi^2];

%L'erreur entre les valeurs propres calculs et analytiques
L=[];
for k=1:size(Vth,2)
L=[L,norm((Vth(k)-Vdisc(k))/Vth(k))];
endfor
% Graphique des fonctions propres calculs
%x=1:10;
%figure 1
%plot(x,V(1),x,V(2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tche 4  6 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cration de la matrice B  partir de A
B=A;

% Numros des lignes/ colonnes  lminer de A
h=[];
for i=N/2:N-1
for j=N/2:N-1

h=[h,i+(N-1)*(j-1)];
endfor
endfor

% Ordonner les numros des lignes/colonnes
h=sort(h);

% Utilisation de la fonction transf qui permet l'limination de la ligne et colonne k
k=0;
for i=1:size(h,2)
j=h(i);
B=transf(B,j-k);
k=k+1;
endfor


% Les 10 plus petites valeurs propres de A
[V1,D1]=eigs(B,10,'sm');
% On met les valeurs prores en ordre croissant
Vdisc1=sort(diag(D1)');
%%les valeurs propores analytiques
Vth=[pi^2 8*pi^2 18*pi^2 2*pi^2 13/4*pi^2 13/4*pi^2 9/2*pi^2 25/4*pi^2 25/4*pi^2 8*pi^2];
%L'erreur entre les valeurs propres calculs et analytiques
L1=[];
for k=1:size(Vth,2)
L1=[L1,norm((Vth(k)-Vdisc1(k))/Vth(k))];
endfor


