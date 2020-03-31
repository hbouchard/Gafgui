
% --------------------------------------------------------------------
function [x,y,nx,ny] = fct_getpoints(z,delta);

[x,y] = ginput;
x=max(min(x,(size(z,2)-1).*delta/2),-(size(z,2)-1).*delta/2);
y=max(min(y,(size(z,1)-1).*delta/2),-(size(z,1)-1).*delta/2);
nx=1+floor((mean(x)+size(z,2).*delta/2)./delta);
ny=1+floor((mean(y)+size(z,1).*delta/2)./delta);