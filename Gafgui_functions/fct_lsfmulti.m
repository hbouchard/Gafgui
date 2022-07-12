% --------------------------------------------------------------------
function [p,sigma0,a1,a2,R,df] = fct_lsfmulti(DOSE,nTHETA,sTHETA,THETA0,M,type)

x = DOSE(:);
y = nTHETA(:);
%Create function matrix and do LSF
F = fct_F_matrix(x,M,type);
[Q,R] = qr(F,0);
%     ws = warning('off','all');
%a = inv(R)*Q'*y;
%a = F\y;
p = inv(F'*F)*F'*y;
r = y - F*p;
df = max(0,length(y) - M);
sigma0 = norm(r)/sqrt(df);

THETA = nTHETA + THETA0;    
p1 = [THETA.^0 THETA.^1]\log(sTHETA);
a1 = p1(1); a2 = p1(2);