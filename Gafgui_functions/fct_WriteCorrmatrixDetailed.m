% --------------------------------------------------------------------
function err = fct_WriteCorrmatrixDetailed(fname,rawres,rawimsize,corrdir,xkrange,rrange,grange,brange,pr,pg,pb);
%the correction works like that
% [XI,SI] = meshgrid(sort(unique(xvect)),[1 2^16-1]);
% [lin,col] = size(XI);
% R0I = zeros(lin,col); G0I = zeros(lin,col); B0I = zeros(lin,col);
% for i=1:col
%     R0I(:,i) = polyval(pr(i,:),SI(:,i));
%     G0I(:,i) = polyval(pg(i,:),SI(:,i));
%     B0I(:,i) = polyval(pb(i,:),SI(:,i));
% end
% r0hat = interp2(XI,SI,R0I,xvect,rvect);
% g0hat = interp2(XI,SI,G0I,xvect,gvect);
% b0hat = interp2(XI,SI,B0I,xvect,bvect);

[~,n] = size(pr);
err = 1;
file = fopen(fname,'w');
fprintf(file,'%.20e\n',rawres);
fprintf(file,'%d\n',rawimsize(1));
fprintf(file,'%d\n',rawimsize(2));
fprintf(file,'%d\n',corrdir);
fprintf(file,'%d\n',min(xkrange));
fprintf(file,'%d\n',max(xkrange));
fprintf(file,'%d\n',min(rrange));
fprintf(file,'%d\n',max(rrange));
fprintf(file,'%d\n',min(grange));
fprintf(file,'%d\n',max(grange));
fprintf(file,'%d\n',min(brange));
fprintf(file,'%d\n',max(brange));
for j=1:n;fprintf(file,'%.20e\n',pr(1,j));end
for j=1:n;fprintf(file,'%.20e\n',pr(2,j));end
for j=1:n;fprintf(file,'%.20e\n',pr(3,j));end
for j=1:n;fprintf(file,'%.20e\n',pg(1,j));end
for j=1:n;fprintf(file,'%.20e\n',pg(2,j));end
for j=1:n;fprintf(file,'%.20e\n',pg(3,j));end
for j=1:n;fprintf(file,'%.20e\n',pb(1,j));end
for j=1:n;fprintf(file,'%.20e\n',pb(2,j));end
for j=1:(n-1);fprintf(file,'%.20e\n',pb(3,j));end
fprintf(file,'%e',pb(3,n))
fclose(file);
err = 0;