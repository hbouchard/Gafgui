function [posx,posy,rmat,bmat,gmat,res,nbbands,stackdir,bits,posxk,posyk,rawimsize] = fct_GetBandsAuto(handles);
%18 March 2020: I added to posxk and posyk, which are indeces of the images
%being croped, being compatible with the detailed homog correction method
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
% button = questdlg('Do you want to use a single film or several films?','Type','Single','Several','Single') ;
button = 'Single';
flag = -1;
if strcmp(button,'Single')
    [ifilename,ipathname]=uigetfile({'*.tif'},'Choose the correction band scan');
    handles.H = fct_updatedisplay(handles);
%     figure(handles.H);
    if ~strcmp(class(ifilename),'double')
        tmpflg  = 1;
    else
        tmpflg  = 0;
    end    
    if tmpflg
        [Im,res,bits] = fct_read_tif_image(fct_makecleanfilename(ipathname,ifilename),'All');
        [lin,col,depth] = size(Im);
        %this is the size of the uncorrected image
        rawimsize = [lin col];
        if depth==1
            Im(:,:,1) = Im;
            Im(:,:,2) = Im(:,:,1);
            Im(:,:,3) = Im(:,:,1);
        elseif depth==2
            Im(:,:,3) = Im(:,:,1);
        end
        figure('NumberTitle','off','Name','Measured bands');
        h0 = fct_display(Im,res);
        h = msgbox('Click on the left-right edges then press enter (only 2 clicks if bands are horizontal)');
        uiwait(h);
        [x,dummy1,nx,dummy2]= fct_getpoints(Im,res);
        if length(x)>1            
            h = msgbox('Click on the up-down edges then press enter (only 2 clicks if bands are vertical)');
            uiwait(h);
            [dummy1,y,dummy2,ny]= fct_getpoints(Im,res);
        else y = [];
        end
        close(h0);
        %if ((length(x)==2)||(length(y)==2))&&((length(x)>2)||(length(y)>2))
        if ((length(x)==2)&&(length(y)>2))||((length(x)>2)&&(length(y)==2))
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
            minxk =  min(min(xgrid))/res;
            minyk =  min(min(ygrid))/res;
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
                %get the indeces
                %18 March 2020: Here I add a patch to get the actual indeces for the
                %detailed homogeneity correction
                for i= 1:nbbands
                    posxk{i} = posx{i}/res - minxk + 1;
                    posyk{i} = posy{i}/res - minyk + 1;
                end
            end
        else
           errordlg('Unable to resolve the edges.'); 
           flag = -1;
        end
    end
else
    tmpflg = 1;
    i = 0;
    restmp = 0;
    bits = 0;
    while tmpflg
        stmp = sprintf('Choose the correction film #%d',i+1);
        [ifilename,ipathname]=uigetfile({'*.tif'},stmp);
        handles.H = fct_updatedisplay(handles);
%         figure(handles.H);
        if ~strcmp(class(ifilename),'double')
            tmpflg  = 1;
        else
            tmpflg  = 0;
        end
        tmpbits = -1;
        if tmpflg
            [Im,res,tmpbits] = fct_read_tif_image(fct_makecleanfilename(ipathname,ifilename),'All');
            figure('NumberTitle','off','Name','Measured bands');
            h = fct_display(Im,res);
            [nlines,ncols] = size(Im(:,:,1));
            %this is the size of the uncorrected image
            [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,res);
            ok = 0;
            if restmp==0
               %this is the first time the variables are defined
               restmp = res;
               rawimsize = [nlines ncols];
               minxk =  min(min(xgrid))/res;
               minyk =  min(min(ygrid))/res;
               ok = 1;
            elseif (res==restmp)&&(rawimsize(1)==nlines)&&(rawimsize(2)==ncols) 
               ok = 1;
            end
            if ok
               [newz,center,owidth,cropz] = fct_getroi(Im,res,'Free',0);
               [m,n,k] = size(newz);
               [GX,GY] = meshgrid(0:n-1,0:m-1);
               i = i+1;
               rmat{i} = newz(:,:,1); gmat{i} = newz(:,:,2); bmat{i} = newz(:,:,3);
               ixstart = round((center(1) - owidth(1)/2)/res);
               iystart = round((center(2) - owidth(2)/2)/res);
               posx{i} = ixstart + GX ;
               posy{i} = iystart + GY ;
            else
                   error('Resolution or image size is not consistent.');
            end
        end
        if bits==0
            bits = tmpbits;
        else
           if  bits~=tmpbits
               error('Bit depth is not consistent between images.');
           end
        end
        button = questdlg('More images?','Shortcut','Yes','No','Yes') ;
        if strcmp(button,'No')
           tmpflg = 0; 
        end
        close ;
    end
    nbbands = i;
    [posx,posy,rmat,bmat,gmat,stackdir] = fct_ReduceBands(posx,posy,rmat,bmat,gmat,nbbands);
    %22 Nov 2016, NPL: here I need to be careful because in the above code
    %I am treating posx as indexes while they should be in fact positions
    %as output of this present function so I fix it by multiplying the
    %indexes by the resolution
    for i=1:nbbands
        posx{i} = posx{i}*res;
        posy{i} = posy{i}*res;
    end
    %18 March 2020: Here I add a patch to get the actual indeces for the
    %detailed homogeneity correction
    for i= 1 :nbbands
        posxk{i} = posx{i}/res - minxk + 1;
        posyk{i} = posy{i}/res - minyk + 1;
    end
    flag = 1;
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
    bits = 0;
end