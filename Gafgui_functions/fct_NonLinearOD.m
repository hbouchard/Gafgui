function [X,m,n] = fct_NonLinearOD(xmat,a,b)

X = [];
for j=1:numel(xmat)
    x = xmat{j} -log10(1+a*10.^(-xmat{j})) -abs(b);
    X = cat(1,X,x);
end
[m,n] = size(X);
X = X(:);