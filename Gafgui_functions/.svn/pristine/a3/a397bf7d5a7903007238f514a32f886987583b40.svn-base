% --------------------------------------------------------------------
function [newz,center,owidth,cropz] = fct_getroi_debug(z,delta,option,iwidth);

zoominopt = 0;
offset = [0 0];
if strcmp(option,'Zoom')
    [z_zoom,rectoff] = imcrop;
    x1 = delta*round(rectoff(1)/delta);
    x2 = delta*round((rectoff(1)+rectoff(3))/delta);
    y1 = delta*round(rectoff(2)/delta);
    y2 = delta*round((rectoff(2)+rectoff(4))/delta);
    offset = [x1 + x2  y1 + y2]/2;
    figure('NumberTitle','off','Name','Zoom');
    hroi = fct_display(z_zoom,delta);
    option = 'Rectangular';
    zoominopt = 1;
elseif strcmp(option,'Zoom-fixed')
    [z_zoom,rectoff] = imcrop;
    x1 = delta*round(rectoff(1)/delta);
    x2 = delta*round((rectoff(1)+rectoff(3))/delta);
    y1 = delta*round(rectoff(2)/delta);
    y2 = delta*round((rectoff(2)+rectoff(4))/delta);
    offset = [x1 + x2  y1 + y2]/2;
    figure('NumberTitle','off','Name','Zoom');
    hroi = fct_display(z_zoom,delta);
    option = 'Fixed';
    zoominopt = 1;
end

if (length(iwidth)<2)&&(~strcmp(option,'Free'))
    width_rect  = inputdlg({'Region width in cm:','Region height in cm:'},'Rectangular region',1,{'1','1'});
    iwidth = [str2double(width_rect(1)) str2double(width_rect(2))];
else
    iwidth = iwidth(:)';
end

if strcmp(option,'Free')
    [temp_newz1,rect] = imcrop;
    x1 = delta*round(rect(1)/delta);
    x2 = delta*round((rect(1)+rect(3))/delta);
    y1 = delta*round(rect(2)/delta);
    y2 = delta*round((rect(2)+rect(4))/delta);
    center = [x1 + x2  y1 + y2]/2;
    center = center + offset;
    owidth = [x2 - x1 y2 - y1];
    %The following doesn't work!!! I doubt there is a weakness here
    [nlines,ncols] = size(z(:,:,1));
    [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,delta);
    xindex = fct_postoindex(center(1)-owidth(1)/2,xgrid);
    yindex = fct_postoindex(center(2)-owidth(2)/2,ygrid);
    cropz = [xindex yindex owidth(1)/delta owidth(2)/delta];
    temp_newz2 = imcrop(z,cropz);
    %it is possible that size(temp_newz1) ~= size(temp_newz2)
    %err = max(abs(size(temp_newz1) - size(temp_newz2)))
elseif  strcmp(option,'Rectangular')||strcmp(option,'Fixed')
    hrect = imrect(gca,[0 0 iwidth(1) iwidth(2)]);
    api = iptgetapi(hrect);
    if strcmp(option,'Fixed')
        api.setResizable(false);
    else
        api.setResizable(true);
    end
    id1 = api.addNewPositionCallback(@(p2) title(mat2str(p2)));
    
%     Changed from waitforbuttonpress to pause. Drag the rectangular then
%     press enter to continue - 3/3/2014 I Billas.
    pause
    
    %     w = 0;
    %     while w==0
    %         w = waitforbuttonpress;
    %     end
    
    position = api.getPosition();
    
    api.removeNewPositionCallback(id1);
    center = [position(1) + position(3)/2 position(2) + position(4)/2];
    center = center + offset;
    owidth = [position(3) position(4)];
    [nlines,ncols] = size(z);
    [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,delta);
    api.delete();
    xindex = fct_postoindex(center(1)-owidth(1)/2,xgrid);
    yindex = fct_postoindex(center(2)-owidth(2)/2,ygrid);
    api.delete();
    cropz = [xindex yindex owidth(1)/delta owidth(2)/delta];
elseif  strcmp(option,'Point')
    [x,y] = ginput;
    center = [mean(x) mean(y)];
    center = [x(length(x)) y(length(y))];
    %ici la fonction n'as pas ete testee avec le zoom, donc la somme du
    %offset de la ligne suivante n'est pas sur de foncitonner dans  imcrop
    %ci-bas
    center = center + offset;
    owidth = iwidth;
    [nlines,ncols] = size(z(:,:,1));
    [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,delta);
    xindex = fct_postoindex(center(1)-owidth(1)/2,xgrid);
    yindex = fct_postoindex(center(2)-owidth(2)/2,ygrid);
    cropz = [xindex yindex owidth(1)/delta owidth(2)/delta];
end
if zoominopt
    close(hroi);
end
newz = imcrop(z,cropz);
imagesc(newz);