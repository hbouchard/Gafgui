% --------------------------------------------------------------------
function [Rot,rrange,grange,brange,trange,err]  = fct_ReadMultiCorrection(fname)

%LUT is a struct with fields rcorr, gcorr and bcorr; 
%the index equals raw signal+1 

err = 1;
Rot = [];
%
file = fopen(fname,'r');
A = fscanf(file,'%f %f %f',[3 inf]);
A = A';
Rot = A(1:3,:);
rrange = A(4,1:2);
grange = A(5,1:2);
brange = A(6,1:2);
trange = A(7,1:2);
[lin,col] = size(Rot);
%
if (size(A,1)==7)&&(size(A,2)==3)
    err = 0;
else
    err = 1
    size(A)
end