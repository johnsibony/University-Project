function C=transf(U,k)

% Enlever les ligne et colonne k de la matrice U et mettre le résultat dans la matrice C

C=[U(1:k-1,1:k-1) U(1:k-1,k+1:end) ; U(k+1:end,1:k-1) U(k+1:end,k+1:end)];
end