% --------------------------------------------------------------------
function LUT  = fct_ReadSignalCorrection(fname)

%LUT is a struct with fields rcorr, gcorr and bcorr; 
%the index equals raw signal+1 

err = 1;
%
file = fopen(fname,'r');
A = fscanf(file,'%d %d %d',[3 inf]);
A = A';
fclose(file);
%this is the uncorrected signal, either red, green or blue
LUT.raw = 0:(2^16-1);
%this is the corrected signal associated to raw
LUT.red = A(:,1);
LUT.green = A(:,2);
LUT.blue = A(:,3);
%
if length(A)==2^16
    err = 0;
else
    err =1
    length(A)
end