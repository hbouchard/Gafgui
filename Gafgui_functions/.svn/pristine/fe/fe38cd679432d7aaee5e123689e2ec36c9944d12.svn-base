% --------------------------------------------------------------------
function H = fct_H_matrix(x,a,type)

N = length(x);
M = length(a);
H = zeros(N,N);

dF = fct_dF_matrix(x,M,type);
h = dF*a;
for i = 1:N
    H(i,i) = h(i);
end