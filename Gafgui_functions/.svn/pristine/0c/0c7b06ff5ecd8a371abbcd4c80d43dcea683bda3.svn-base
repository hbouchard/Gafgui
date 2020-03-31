% --------------------------------------------------------------------
function [var,profil] = fct_getmeanprofilefrompoint(z,delta,dir,shift)


[x,y,nx,ny] = fct_getpoints(z,delta);
if strcmp(dir,'x')==1
    k1 = max((ny-shift),1);
    k2 = min((ny+shift),size(z,1));
    profil = mean(z(k1:k2,:),1);
    var=(1:size(z,2))*delta;
elseif strcmp(dir,'y')==1
    k1 = max((ny-shift),1);
    k2 = min((ny+shift),size(z,2));
    profil = mean(z(:,k1:k2),2);
    var=(1:size(z,1))*delta;
else
    error('Wrong direction option.');
end
