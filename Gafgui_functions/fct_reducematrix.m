% --------------------------------------------------------------------
function Y = fct_reducematrix(X,f);

[M,N] = size(X);
M = floor(M/f)*f;
N = floor(N/f)*f;
m = M/f;
n = N/f;
Y = zeros(m,n);
for i=1:f
    for j=1:f
        Y = Y + X((0:(m-1))*f+i,(0:(n-1))*f+j);
    end
end
Y = Y/f^2;