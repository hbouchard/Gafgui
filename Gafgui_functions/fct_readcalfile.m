% --------------------------------------------------------------------
function [DOSE,OD,M,type,sigparam,Npix,THETA0] = fct_readcalfile(file)

%WB july 2022: added bck values, Nbackground value and s0 (see end of fct_CreateCalCurveNewMultiMethod)
A = fscanf(file, '%e\t%e\t%e\n', [3, inf]);
A = A';
k = size(A,1);
DOSE = A(1:k-3,1);
DOSE = abs(DOSE);
OD = A(1:k-3,2);
s1 = A(1:k-3,3);
s0 = A(k-2,3);
sigparam{1} = s0; sigparam{2} = s1;
Npix(1) = A(k-2, 1);
Npix(2) = A(k-2, 2);
type = A(k-1,3);
M = A(k-1,2);
THETA0 = A(k,1);

end



% 
% Npix(1:2) = 0;
% sigparam(1:2) = 0;
% A = fscanf(file,'%f\t%f',[2 inf]);
% A = A';
% k = size(A,1);
% DOSE = A(:,1);
% 
% 
% DOSE = abs(DOSE);
% OD = A(:,2);
% type = DOSE(k);
% M = OD(k);
% k = k-1;
% DOSE = DOSE(1:k);
% OD = OD(1:k);
% 
% sigparam(1) = DOSE(k);
% sigparam(2) = OD(k);
% k = k-1;
% DOSE = DOSE(1:k);
% OD = OD(1:k);
% 
% Npix(1) = DOSE(k);
% Npix(2) = OD(k);
% k = k-1;
% DOSE = DOSE(1:k);
% OD = OD(1:k);