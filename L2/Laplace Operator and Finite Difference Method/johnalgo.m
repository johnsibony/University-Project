
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%˘
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T‚che 1 ‡ 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Nombre de points dans la grille
N=10;
% Pas de la grille
h=2/N;


% Les quatres listes de numÈrotations Nord/ Sud/Est et Ouest
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NumÈros des lignes/ colonnes ‡ ÈlÈminer de A
a=[];
for i=N/2:N-1
for j=N/2:N-1

a=[a,i+(N-1)*(j-1)];
endfor
endfor

% Ordonner les numÈros des lignes/colonnes
a=sort(a);

% Liste Nord2 et Sud2
Nord2=[N:(N-1)*(N-1) zeros(1,N-1)];
Sud2=[zeros(1,N-1) 1:(N-1)*(N-2)];
for k=1: (N-1)*(N-1)
for j=1:N*N/4
if Nord2(k)==a(j)
Nord2(k)=0;
endif
if Sud2(k)==a(j)
Sud2(k)=0;
endif
endfor
endfor

%Liste Ouest2
Ouest2=[0 1:(N-1)*(N-1)-1];
for k=1:N-2
j=(N-1)*k+1;
Ouest2(j)=0;
endfor
for k=1: (N-1)*(N-1)
for j=1:N*N/4
if Ouest2(k)==a(j)
Ouest2(k)=0;
endif
endfor
endfor

%Liste Est2
Est2=[2:(N-1)*(N-1) 0];
for k=1:N-2
j=(N-1)*k;
Est2(j)=0;
endfor
for k=1: (N-1)*(N-1)
for j=1:N*N/4
if Est2(k)==a(j)
Est2(k)=0;
endif
endfor
endfor

% Remplissage de la diagonale
B=diag(4/h^2*ones(1,(N-1)*(N-1)));

% Format creuse
B=sparse(B);

% Utilisation de la liste Nord 
for k=1:(N-1)*(N-1)
if Nord2(k)~=0
B(k,Nord2(k))=-1/h^2;
endif

%Utilisation de la liste Sud
if Sud2(k)~=0
B(k,Sud2(k))=-1/h^2;
endif

%Utilisation de la liste Est
if Est2(k)~=0
B(k,Est2(k))=-1/h^2;
endif

%Utilisation de la liste Ouest
if Ouest2(k)~=0
B(k,Ouest2(k))=-1/h^2;
endif

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

