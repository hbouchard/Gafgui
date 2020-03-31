function r = fct_minrealcubicroot(c3,c2,c1,c0)
% 
% clc;
% clear all;
% close all;
% 
% c3 = 1; c2 = -2; c1 = -1; c0 = 2;

r = fct_real_cubic_roots(c3,c2,c1,c0);
%find the nonreal (or zero) roots
k1 = find(r(:,1)==0);
k2 = find(r(:,2)==0);
k3 = find(r(:,3)==0);
%put a very large number in these ones
r(k1,1) = r(k1,1) +1e10;
r(k2,2) = r(k2,2) +2e10;
r(k3,3) = r(k3,3) +3e10;

%then test with large numbers instead of zeros
%find the minimum
%I is the column index where the min is
[dummy,j] = min(abs(r),[],2);
[m,n] = size(r); %here the index is i,j
r = reshape(r',m*n,1); %here the index is k=n*(i-1)+j
r = r(n*((1:m)'-1)+j);

%put back zeros in large numbers. they only exists if there was no smaller
%ones
%if there is no nonreal root, just put it to 0
r(find(r>9e9)) = 0;



