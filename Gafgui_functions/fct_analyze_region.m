% --------------------------------------------------------------------
function [ave, stdev, N] = fct_analyze_region(A)

% ave = mean(mean(A));
% B = (A-ave).^2;
% A = A.*0+1;
% N = sum(sum(A));
% stdev = sum(sum(B))/(N-1);
% stdev = sqrt(stdev);

[m,n] = size(A);

ave = mean(reshape(A,m*n,1));
stdev = std(reshape(A,m*n,1));
N = m*n;