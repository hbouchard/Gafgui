% --------------------------------------------------------------------
function [odi,dosei] = fct_getdose_from_od(DOSE,OD,M,type)

[p,sy,R,df] = fct_lsf(DOSE,OD,M,type);
dosei = 0:1:max(DOSE);
dosei = dosei(:);
Fi = fct_F_matrix(dosei,M,type);
odi = Fi*p;
Q = Fi*inv(R);
%maniere pratique de resoudre puisque la metode qr(F,0)
%done une matrice R carree, on peut resoudre V = F*inv(R'*R)*F'
%par V  = F*inv(R)*(F*inv(R))'=Q*Q'. Les elements diag sont
%donc la lgine suivante
sodi = sy*sqrt(sum(Q.*Q,2));
%%%%%%
dFi = fct_dF_matrix(dosei,M,type);
d2Fi = fct_d2F_matrix(dosei,M,type);
bias = d2Fi*p/2./(dFi*p).^3.*(sodi.^2-sy^2);
%%%%%%
dosei = dosei(:)-bias(:);