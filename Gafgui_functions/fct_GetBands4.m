function [posx,posy,rmat,bmat,gmat,res,nbbands,stackdir] = fct_GetBands4(Im,res);
% % 
% clc;
% clear all;
% close all;
% 
% % This add the functions path
% cpu = computer;
% if strcmp(cpu,'PCWIN')||strcmp(cpu,'PCWIN64')
%     c = '\'; % dos
% else
%     c = '/'; % mac/linux/unix
%     %In the case someone uses a SUN, I don't know if this is correct
% end
% path = cd;
% if (path(max(size(path,1),size(path,2)))==c)
%     path = [path 'Gafgui_functions'];
%     path = [path '..' c 'Gafgui']
% else
%     path = [path c 'Gafgui_functions'];
%     path = [path c '..' c 'Gafgui'];
% end
% addpath(path);
% fct_AddGafguiFctPath();
% handles = fct_initGafgui();

%%%%%%%%%%%%%%%%%%%%%%%
%[Im,res] = fct_read_tif_image(fct_makecleanfilename(ipathname,ifilename),'All');
figure('NumberTitle','off','Name','Correction bands');
h0 = fct_display(Im,res);
h = msgbox('Click on all the horizontal edges then press enter');
uiwait(h);
[x,dummy1,nx,dummy2]= fct_getpoints(Im,res);
if length(x)>1            
    h = msgbox('Click on all the vertical edges then press enter');
    uiwait(h);
    [dummy1,y,dummy2,ny]= fct_getpoints(Im,res);
else y = [];
end
close(h0);
if ((length(x)==2)||(length(y)==2))&&((length(x)>2)||(length(y)>2))
    x = sort(x); y = sort(y);
    nbbands = max(length(x),length(y))-1;
    if length(y)>length(x)
       stackdir = 'Vertical';
    else
       stackdir = 'Horizontal';  
    end
    for i= 1 :nbbands 
        if strcmp(stackdir,'Vertical')
            xdelim{i} = [x(1) x(2)];
            ydelim{i} = [y(i) y(i+1)];
        elseif strcmp(stackdir,'Horizontal')
            ydelim{i} = [y(1) y(2)];
            xdelim{i} = [x(i) x(i+1)];
        end
    end
    [nlines,ncols] = size(Im(:,:,1));
    [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,res);
    xgrid = xgrid(:);
    ygrid = ygrid(:);
    flag = 1;
    while flag==1
        d = inputdlg({'Enter the distance in mm from edges to eliminate'},'Eliminate',1,{'5'}) ;
        if length(d)==0
            flag = -1;
        else
            d = str2double(d);
            d = max(0,d)/10;
            im = [];
            x = [];
            y = [];
            for i= 1 :nbbands
               kx1 = find(abs(xdelim{i}(1)+d-xgrid)==min(abs(xdelim{i}(1)+d-xgrid)));
               kx2 = find(abs(xdelim{i}(2)-d-xgrid)==min(abs(xdelim{i}(2)-d-xgrid)));
               ky1 = find(abs(ydelim{i}(1)+d-ygrid)==min(abs(ydelim{i}(1)+d-ygrid)));
               ky2 = find(abs(ydelim{i}(2)-d-ygrid)==min(abs(ydelim{i}(2)-d-ygrid)));
               if strcmp(stackdir,'Vertical')
                    im = cat(1,im,Im(ky1:ky2,kx1:kx2,:));
               elseif strcmp(stackdir,'Horizontal')
                    im = cat(2,im,Im(ky1:ky2,kx1:kx2,:));                   
               end
               x = cat(1,x,xgrid(kx1:kx2));
               y = cat(1,y,ygrid(ky1:ky2));                       
            end
            h = figure;
            imagesc(x,y,im);
            set(gca,'DataAspectRatio',[1 1 1]);
            set(gcf,'Units','pixels');
            button = questdlg('Do you want to keep these bands?','Bands','Yes','No','Yes') ;
            if strcmp(button,'Yes')
                flag = 0;
            elseif strcmp(button,'No')
                flag = 1;
            else
                flag = -1;
            end
            close(h);
        end
    end
    if flag ~=-1
        im = [];
        x = [];
        y = [];
        for i= 1 :nbbands
           kx1 = find(abs(xdelim{i}(1)+d-xgrid)==min(abs(xdelim{i}(1)+d-xgrid)));
           kx2 = find(abs(xdelim{i}(2)-d-xgrid)==min(abs(xdelim{i}(2)-d-xgrid)));
           ky1 = find(abs(ydelim{i}(1)+d-ygrid)==min(abs(ydelim{i}(1)+d-ygrid)));
           ky2 = find(abs(ydelim{i}(2)-d-ygrid)==min(abs(ydelim{i}(2)-d-ygrid)));
           rmat{i} = Im(ky1:ky2,kx1:kx2,1);
           gmat{i} = Im(ky1:ky2,kx1:kx2,2);
           bmat{i} = Im(ky1:ky2,kx1:kx2,3);
           [posx{i},posy{i}] = meshgrid(xgrid(kx1:kx2),ygrid(ky1:ky2)) ;
           if strcmp(stackdir,'Vertical')
                im = cat(1,im,Im(ky1:ky2,kx1:kx2,:));
           elseif strcmp(stackdir,'Horizontal')
                im = cat(2,im,Im(ky1:ky2,kx1:kx2,:));                   
           end
           x = cat(1,x,xgrid(kx1:kx2));
           y = cat(1,y,ygrid(ky1:ky2));                        
        end
%                 h = figure;
%                 imagesc(x,y,im);
%                 set(gca,'DataAspectRatio',[1 1 1]);
%                 set(gcf,'Units','pixels');
    end
else
   errordlg('Unable to resolve the edges.'); 
   flag = -1;
end

if flag==-1
    posx = [];
    posy = [];
    rmat = [];
    gmat = [];
    bmat = [];
    res = 0;
    nbbands = 0;
    stackdir = [];
end