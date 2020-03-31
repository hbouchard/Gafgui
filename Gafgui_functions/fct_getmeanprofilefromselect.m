% --------------------------------------------------------------------
function [VAR,PROFIL,newz] = fct_getmeanprofilefromselect(z,delta,dir,button,channel,iwidth);

[newz,center,width] = fct_getroi(z(:,:,channel),delta,button,iwidth);
%newz = newz(:,:,ichannel);
[m,n] = size(newz);
VAR = [];
PROFIL = [];
if m*n~=0
    for ichannel=channel
        if strcmp(dir,'x')==1
            profil = mean(newz(:,:,ichannel),1);
            var = ((0:size(profil,2)-1)-0.5*(size(profil,2)-1))*delta;
            var = var + center(1);
        elseif strcmp(dir,'y')==1
            profil = mean(newz(:,:,ichannel),2);
            var=((0:size(profil,1)-1)-0.5*(size(profil,1)-1))*delta;
            var = var + center(2);
        end    
        VAR = cat(2,VAR,var(:));
        PROFIL= cat(2,PROFIL,profil(:));
    end
end