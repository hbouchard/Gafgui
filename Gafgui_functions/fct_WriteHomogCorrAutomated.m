% --------------------------------------------------------------------
function err = fct_WriteHomogCorrAutomated(fname,PR,PG,PB,x,rlims,glims,blims)

nR = length(PR(:));
nG = length(PG(:));
nB = length(PB(:));
err = 1;
file = fopen(fname,'w');
fprintf(file,'%d\n',rlims(1));
fprintf(file,'%d\n',rlims(2));
fprintf(file,'%d\n',glims(1));
fprintf(file,'%d\n',glims(2));
fprintf(file,'%d\n',blims(1));
fprintf(file,'%d\n',blims(2));
fprintf(file,'%d\n',nR);
fprintf(file,'%d\n',nG);
fprintf(file,'%d\n',nB);
for j=1:nR;fprintf(file,'%e\n',PR(j));end
for j=1:nG;fprintf(file,'%e\n',PG(j));end
for j=1:nB;fprintf(file,'%e\n',PB(j));end
for j=1:(length(x)-1);fprintf(file,'%f\n',x(j));end
fprintf(file,'%f',x(end));

fclose(file);
err = 0;