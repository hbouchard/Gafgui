% --------------------------------------------------------------------
function [corrres,rawimsize,corrdir,xkrange,rrange,grange,brange,pr,pg,pb] = fct_ReadCorrmatrixDetailed(fname)

%the correction works like that
% [XI,SI] = meshgrid(unique(xvect),[1 2^16-1]);
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

file = fopen(fname,'r');
A = fscanf(file,'%e',[1 inf]);
fclose(file);

A = A(:);
corrres = A(1);
corrres = 2.54/round(2.54/corrres);
%Here I need to put indeces rather than positions. The reason is round off
%errors occuring while you read the input file. The variable res is
%typically 2.54/DPI, where DPI is an integer.
rawimsize = [A(2) A(3)];
corrdir = A(4);
xkrange = [A(5) A(6)];
n = xkrange(2)-xkrange(1)+1;
rrange = [A(7) A(8)];
grange = [A(9) A(10)];
brange = [A(11) A(12)];
A = A(13:end);
if (length(A)/9)~=n%mod(length(A),9)~=0
    error('Wrong data length, should be an integer factor of 9 for quadratic fits.');
else
    pr = [A(1:n) A(n+1:2*n) A(2*n+1:3*n)]';
    A = A(3*n+1:end);
    pg = [A(1:n) A(n+1:2*n) A(2*n+1:3*n)]';
    A = A(3*n+1:end);
    pb = [A(1:n) A(n+1:2*n) A(2*n+1:3*n)]';
end

err = 0;