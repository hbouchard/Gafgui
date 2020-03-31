% --------------------------------------------------------------------
function i = fct_postoindex(x,xgrid)

N = max(size(xgrid));
delta = abs(xgrid(1)-xgrid(2));
i = (x - xgrid(1))/delta +1 ;
i = min(max(round(i),1),N);