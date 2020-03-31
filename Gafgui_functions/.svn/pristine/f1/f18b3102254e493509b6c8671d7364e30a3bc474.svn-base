function [a,b] = fct_FindValuesMinRowColMatrices(A,B,dim)
% 
% clc;
% clear all;
% close all;
% format short
% A  =rand(3,2)
% B  =[1 2; 3 4; 5 6;]
% dim = 1;

[dummy,jmin] = min(A,[],dim); 
%here the index is i,j
[lin,col] = size(A);
if dim ==2
    a = reshape(A',lin*col,1); %here the index is k=col*(i-1)+j
    b = reshape(B',lin*col,1);
    kmin = col*((1:lin)'-1)+jmin; 
elseif dim==1
    a = reshape(A,1,lin*col); %here the index is k=n*(i-1)+j
    b = reshape(B,1,lin*col);
    kmin = lin*((1:col)-1)+jmin;
end

a = a(kmin);
b = b(kmin);
