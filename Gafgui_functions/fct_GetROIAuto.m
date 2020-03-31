function [posx,posy,rmat,bmat,gmat,res,nbROI,bits] = fct_GetROIAuto(handles);
% % 
% clc;
% clear all;
% close all; 
% fct_AddGafguiFctPath();
% handles = fct_initGafgui();

flag =1;
[ifilename1,ipathname1]=uigetfile({'*.tif'},'Choose unirradiated films scan');
if (strcmp(class(ifilename1),'double'))
    flag = 0;
else    
    [ICAL,res,bits] = fct_read_tif_image(fct_makecleanfilename(ipathname1,ifilename1),'All');
    figure('NumberTitle','off','Name','Calibration curve: background');
    h = fct_display(ICAL,res);
    ans = questdlg('Do you want to zoom to a smaller ROI?','Image','Yes','No','Yes') ;
    if strcmp(ans,'Yes')
        figure(h);
        ICAL = imcrop;
        h = fct_display(ICAL,res);
    end
    figure(h);
    nbfilms  = inputdlg({'Number of ROI'},'Number of films',1);
    if numel(nbfilms)~=0
        nbfilms  = str2double(nbfilms);
        nbROI = nbfilms;
%         Im = -log10(max(double(ICAL(:,:,1)),1)./(2^(bits)-1));
%         Im = cat(3,Im,-log10(max(double(ICAL(:,:,2)),1)./(2^(bits)-1)));
%         Im = cat(3,Im,-log10(max(double(ICAL(:,:,3)),1)./(2^(bits)-1)));
        [nlines,ncols] = size(ICAL(:,:,1));
        [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,res);
        xgrid = xgrid(:);
        ygrid = ygrid(:);
        type = 'Fixed';
        rectsize  = inputdlg({'Region width in cm:','Region height in cm:'},'Background',1,{'1','1'}) ;
        if numel(rectsize)==0
            flag = 0;
        else
            width_bck = str2double(rectsize(1));
            height_bck = str2double(rectsize(2));
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            BCKGRND = zeros(nbfilms,3);
            NBCK = width_bck*height_bck/(handles.CCDres^2);
            figure(h);
            for i=1:nbfilms
                [newz,center,owidth,cropz] = fct_getroi(ICAL(:,:,1),res,type,[width_bck height_bck]);
                rmat{i} = newz;
                newz = imcrop(ICAL(:,:,2),cropz);
                gmat{i} = newz;
                newz = imcrop(ICAL(:,:,3),cropz);
                bmat{i} = newz;
                [nlines,ncols] = size(newz);
                [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,res);
                [posx{i},posy{i}] = meshgrid(xgrid,ygrid);
            end
            figure(h);
            close(h);
            %%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end
if flag==0
    posx = [];
    posy = [];
    rmat = [];
    gmat = [];
    bmat = [];
    res = 0;
    nbROI = 0;
    bits = 0;
end