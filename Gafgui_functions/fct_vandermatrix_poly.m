% --------------------------------------------------------------------
function [V,nkernels] = fct_vandermatrix_poly(x,order,zero)
V = [];
if zero
    V = 1+x(:)*0;
end
for n = 1:order
    V = cat(2,V,x(:).^n);
end
nkernels = size(V,2);