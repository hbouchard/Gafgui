function [DOSE,nTHETA,sTHETA,THETA0,Npix,res,channel,opt,N] = fct_ReadCalFileMulti(file);

% clear all;
% clear functions;
% close all;
% clc;
% fct_AddGafguiFctPath();
% handles = fct_initGafgui('Black');
% 
% [ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');
% file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');

% %HB: this is Bouchard et al 2021
% if ~strcmp(class(ofilename),'double')
%     file = fopen(fct_makecleanfilename(opathname,ofilename),'w');
%     for i=1:nbfilms
%         fprintf(file,'%e\t%e\t%e\n',DOSE(i),OD(i),sOD(i));
%     end
%     %HB 18 jan 2021: we no longer separately consider the ucnertainty on the background value, i.e. it is built in sig0
%     fprintf(file,'%e\t%e\t%e\n',p1(1), p1(2),THETA0);
%     fprintf(file,'%e\t%e\t%e\n',s0,RES_RAD,Npix);
%     fprintf(file,'%e\t%e\t%e',opt,N,Channel);
%     fclose(file);
% end
            
A = fscanf(file,'%e\t%e\t%e',[3 inf]);
A = A';

channel = A(end,1);
N = A(end,2);
opt = A(end,3);
A = A(1:end-1,:);

THETA0 = A(end,1);
Npix = A(end,2);
res = A(end,3);
A = A(1:end-1,:);

DOSE = A(:,1);
nTHETA = A(:,2);
sTHETA = A(:,3);
