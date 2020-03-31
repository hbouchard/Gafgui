function [DOSE,XI,sXI,c,V,rescal,sximesh] = fct_ReadCalFileMulti(file);

% clear all;
% clear functions;
% close all;
% clc;
% fct_AddGafguiFctPath();
% handles = fct_initGafgui('Black');

%         file = fopen(fct_makecleanfilename(opathname,ofilename),'w');
%         fprintf(file,'%.10e\t%.10e\t%.10e\n',c(1),c(2),RES_RAD);
%         fprintf(file,'%.10e\t%.10e\t%.10e\n',V(1,1),V(1,2),V(2,2));       
%         fprintf(file,'%.10e\t%.10e\t%.10e\n',m,n,sximesh.npix0);
%         for i=1:m
%             for j=1:n
%                 fprintf(file,'%.10e\t%.10e\t%.10e\n',sximesh.x(i,j),sximesh.y(i,j),sximesh.z(i,j));
%                 
%             end
%         end
%         for i=1:length(DOSE)-1
%             fprintf(file,'%.10e\t%.10e\t%.10e\n',DOSE(i),XI(i),sXIhat0(i));
%         end
%         fprintf(file,'%.10e\t%.10e\t%.10e',DOSE(i+1),XI(i+1),sXIhat0(i+1));
          
% [ifilename,ipathname] = uigetfile({'*.mlt'},'Choose calibration curve');
% file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');

A = fscanf(file,'%e\t%e\t%e',[3 inf]);
A = A';
%Entry 1
a = A(1,:); 
c = a(1:2)'; rescal = a(3);
%Entry 2
a = A(2,:);
V = zeros(2,2);
V(1,1) = a(1); V(2,1) = a(2); V(1,2) = a(2); V(2,2) = a(3); 
%Entry 3
a = A(3,:);
m = a(1); n = a(2); npix0 = a(3);
A = A(4:end,:);
%Entry 4
X = []; Y = []; Z = [];
for j=1:n
    a = A(1:m,:);
    X = cat(2,X,a(:,1)); Y = cat(2,Y,a(:,2)); Z = cat(2,Z,a(:,3));
    A = A((m+1):end,:);
end
%Entry 5
a = A(:,:);
DOSE = a(:,1); XI = a(:,2); sXI = a(:,3);

sximesh.npix0 = npix0; sximesh.x = X; sximesh.y = Y; sximesh.z = Z;