function CV = fct_PCACovarMatrix(V)

[N,M] = size(V);
uN = ones(N,1);
uM = ones(M,1);
CV = 1/(N-1)*V'*V - 1/N/(N-1)*V'*uN*uN'*V;