% --------------------------------------------------------------------
function [odi,dosei] = fct_getcalcurvepoints(DOSE,OD,M,type)

DOSE = abs(DOSE);
[p,sy,R,df] = fct_lsf(DOSE,OD,M,type);
dosei = 0:1:floor(max(DOSE));
dosei = dosei(:);
Fi = fct_F_matrix(dosei,M,type);
odi = Fi*p;