function CX = fct_CovarMatrix(X)

[N,M] = size(X);
uN = ones(N,1);
uM = ones(M,1);
CX = 1/(N-1)*X'*X - 1/N/(N-1)*X'*uN*uN'*X;