% --------------------------------------------------------------------
function [I,Sr,Sg,Sb] = fct_ReadODfile(file)

A = fscanf(file,'%f\t%f\t%f\t%f',[4 inf]);
A = A';
k = size(A,1);
I = A(:,1);
Sr = A(:,2);
Sg = A(:,3);
Sb = A(:,4);