

% --------------------------------------------------------------------
function [k,y] = fct_indexreducerand(x,N)

K = min(N,length(x));
if K==N
    rdn = rand(K,1);
    [rdn,k] = sort(rdn);
    k = unique(round(length(x)/K*k));
    k = k(find(k~=0));
    y = x(k);
else
    k = 1:length(x);
    y = x;
end