% --------------------------------------------------------------------
function [ave, stdev, N] = fct_analyze_region_circular(A)

[lin,col] = size(A);
lin0 = (lin-1)/2;
col0 = (col-1)/2;
[j,i] = meshgrid(1:col,1:lin);
i = (i-1-lin0)/min(lin0,col0);
j = (j-1-col0)/min(lin0,col0);
Mask = sqrt(i.^2+j.^2);
Mask = logical(1.0000001-min(Mask,1.0000001));
A = Mask.*A;
N = sum(sum(Mask));
ave = sum(sum(A))/N;
B = (A-ave.*Mask).^2;
stdev = sum(sum(B))/(N-1);
stdev = sqrt(stdev);