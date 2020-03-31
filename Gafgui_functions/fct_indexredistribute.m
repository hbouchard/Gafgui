
% --------------------------------------------------------------------
function [k,y] = fct_indexredistribute(x,N)

bins = (0:N-1)*(max(x)-min(x))/(N-1) + min(x);
k = 1;
for i=1:N-1
    kk = intersect(find(x>bins(i)),find(x<bins(i+1)));
    if length(kk)~=0
        [dmb,ii] = sort(rand(length(kk),1));
        k = cat(1,k,kk(ii(1)));
    end
end
k = cat(1,k,length(x));
y = x(k);