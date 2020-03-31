% --------------------------------------------------------------------
function [a,sigma,R,df] = fct_lsf(DOSE,OD,M,type)

x = DOSE(:);
y = OD(:);
%Create function matrix and do LSF
F = fct_F_matrix(x,M,type);
[Q,R] = qr(F,0);
%     ws = warning('off','all');
%a = inv(R)*Q'*y;
%a = F\y;
a = inv(F'*F)*F'*y;
r = y - F*a;
df = max(0,length(y) - M);
sigma = norm(r)/sqrt(df);
