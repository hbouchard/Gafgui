% --------------------------------------------------------------------
function [DOSE,OD,M,type,sigparam,Npix] = fct_readcalfile(file)

Npix(1:2) = 0;
sigparam(1:2) = 0;
A = fscanf(file,'%f\t%f',[2 inf]);
A = A';
k = size(A,1);
DOSE = A(:,1);


DOSE = abs(DOSE);
OD = A(:,2);
type = DOSE(k);
M = OD(k);
k = k-1;
DOSE = DOSE(1:k);
OD = OD(1:k);

sigparam(1) = DOSE(k);
sigparam(2) = OD(k);
k = k-1;
DOSE = DOSE(1:k);
OD = OD(1:k);

Npix(1) = DOSE(k);
Npix(2) = OD(k);
k = k-1;
DOSE = DOSE(1:k);
OD = OD(1:k);