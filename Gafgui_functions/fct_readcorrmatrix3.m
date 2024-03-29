% --------------------------------------------------------------------
function [p,xmin,xmax,x0,smin,smax] = fct_readcorrmatrix3(fname)

file = fopen(fname,'r');
A = fscanf(file,'%e',[1 inf]);
fclose(file);
type = A(1);
A = A(2:end);
xmin = A(1);
xmax = A(2);
x0 = A(3);
smin = A(4);
smax = A(5);
npol = A(6);
A = A(7:length(A));
orders = A(1:npol);
A = A((npol+1):length(A));
A = A(:);
for j=1:npol
    l = orders(j)/2+1;
    p{j} = A(1:l);
    if j<npol
        A = A((l+1):length(A));
    end
end