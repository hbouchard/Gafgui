% --------------------------------------------------------------------
function err = fct_WriteSignalCorrection(fname,rcorr,gcorr,bcorr)

err = 1;
%
file = fopen(fname,'w');
for i=1:(2^16-1)
    fprintf(file,'%d %d %d\n',uint16(rcorr(i)),uint16(gcorr(i)),uint16(bcorr(i)));
end
fprintf(file,'%d %d %d',uint16(rcorr(end)),uint16(gcorr(end)),uint16(bcorr(end)));
fclose(file);
%
err = 0;