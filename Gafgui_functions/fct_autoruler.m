
% --------------------------------------------------------------------
function res = fct_autoruler(A,delta)

lengthy = size(A,1)*delta;
lengthx = size(A,2)*delta;
res = (lengthx+lengthy)/2/10;
N = ceil(log10(res))-1;
res = 10^N * round(res/10^N);