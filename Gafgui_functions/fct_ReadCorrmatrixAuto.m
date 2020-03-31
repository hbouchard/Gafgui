% --------------------------------------------------------------------
function [xx,s1,s2,s3,c1,c2,c3,p1,p2,p3,order1,order2,sym,m,n,bits] = fct_WriteCorrmatrixAuto(fname)

file = fopen(fname,'r');
A = fscanf(file,'%e',[1 inf]);
fclose(file);
A = A(:);
bits = A(1);
order1 = A(2);
order2 = A(3);
sym = A(4);
m = A(5);
n = A(6);
l = A(7); A = A(8:end);
p1 = A(1:l); A = A(l+1:end);
p2 = A(1:l); A = A(l+1:end);
p3 = A(1:l); A = A(l+1:end);
l = A(1); A = A(2:end);
xx = A(1:l); A = A(l+1:end);
s1 = A(1:l); A = A(l+1:end);
s2 = A(1:l); A = A(l+1:end);
s3 = A(1:l); A = A(l+1:end);
c1 = A(1:l); A = A(l+1:end);
c2 = A(1:l); A = A(l+1:end);
c3 = A(1:l); 
err = 0;