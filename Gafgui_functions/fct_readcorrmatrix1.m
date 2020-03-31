% --------------------------------------------------------------------
function [p,center,xmin,xmax,smin,smax] = fct_readcorrmatrix1(fname)

file = fopen(fname,'r');
A = fscanf(file,'%e',[1 inf]);
fclose(file);
type = A(1);
A = A(2:end);
npol = A(1);
orders = A(2:(npol+1));
xmin = A(npol+2);
xmax = A(npol+3);
smin = A(npol+4);
smax = A(npol+5);
A = A(npol+6:length(A));
A = A(:);
center(1:npol) = 0;
index(1:npol) = 0;
index(1) = 1;
for i = 2:npol
    index(i) = index(i-1) + orders(i-1)+1;
end
for i = 1:(npol-1)
    p{i} = A(index(i):(index(i+1)-1));
end
p{npol} = A(index(npol):length(A));
for i = 1:npol
    center(i) = polyval(p{i},0);
end
center = center(:);