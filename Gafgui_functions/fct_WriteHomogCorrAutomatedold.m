% --------------------------------------------------------------------
function err = fct_WriteHomogCorrAutomated(fname,IXr,ISr,IR0,IXb,ISb,IG0,IXg,ISg,IB0)

%the correction works like that
%     Imcorr(:,:,1) = uint16(interp2(XIr,SIr,R0I,X,R));
%     Imcorr(:,:,2) = uint16(interp2(XIg,SIg,G0I,X,G));
%     Imcorr(:,:,3) = uint16(interp2(XIb,SIb,B0I,X,B));

nxr = length(unique(IXr(:)));
nxg = length(unique(IXg(:)));
nxb = length(unique(IXb(:)));

nsr = length(unique(ISr(:)));
nsg = length(unique(ISg(:)));
nsb = length(unique(ISb(:)));


err = 1;
file = fopen(fname,'w');
fprintf(file,'%d %d %d\n',res);

[lin,col] = size(XI);

for j=1:n;fprintf(file,'%e\n',ix(j));  end
for j=1:n;fprintf(file,'%e\n',pr(j,1));end
for j=1:n;fprintf(file,'%e\n',pr(j,2));end
for j=1:n;fprintf(file,'%e\n',pg(j,1));end
for j=1:n;fprintf(file,'%e\n',pg(j,2));end
for j=1:n;fprintf(file,'%e\n',pb(j,1));end
for j=1:(n-1);fprintf(file,'%e\n',pb(j,2));end
fprintf(file,'%e',pb(n,2))
fclose(file);
err = 0;