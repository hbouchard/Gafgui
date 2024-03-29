
% --------------------------------------------------------------------
function [N,xi,F,center] = fct_readcorrmatrix2(fname)

file = fopen(fname,'r');
A = fscanf(file,'%f',[1 inf]);
fclose(file);
type = A(1);
A = A(2:end);
N = A(1); %number of Y arrays
A = A(2:length(A));
M = length(A); %total number of data
n = M/(N+1); %number of values in x
xi = A(1:n);%array x
A = A((n+1):length(A));
for i=1:N
    F{i} = A(1:n);
    A = A((n+1):length(A));
end
center(1:N) = 0;
for i=1:N
    y = F{i};
    k = min(find(abs(xi)==min(abs(xi))));
    center(i) = y(k);
end