% --------------------------------------------------------------------
function [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,delta)

xgrid = (-(ncols-1)/2 :(ncols-1)/2 )*delta;
ygrid = (-(nlines-1)/2:(nlines-1)/2)*delta;