function [Y,P,CY] = fct_GetPCA(X,normal)

[N,M] = size(X);
%unitary vectors
CX = fct_PCACovarMatrix(X);
%CW = Vw;
%eigenvalues
[E,eigval] = eig(CX);
%I cross my fingers that all eigenvalues are positive... 
P = E';
%%%%%%%%%%%%%%%%%%%%
%normalization: causes a lot of problems
if normal %1 is normalize=yes
    c = 1./sum(P,2);
    C = [];
    for i=1:M
        C = cat(2,C,c);
    end
    P = P.*C;
else %you choose not to normalize
    P = E';
end
%%%%%%%%%%%%%%%%%%%%
%We need to find the right orders of the PC
Y = X*inv(P);
CY = fct_PCACovarMatrix(Y);
[~,isort] = sort(diag(CY));
isort = flipud(isort(:));
Pnew = [];
for i=1:M
   p = P(isort(i),:);
   Pnew = cat(1,Pnew,p);
end
P = Pnew; 
Y = X*inv(P);