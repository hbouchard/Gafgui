% function err = fct_RegisterImage(handles);

clear all;
clear functions;
close all;
clc;
fct_AddGafguiFctPath();
handles = fct_initGafgui();

err = 0;

[ifilename1,ipathname1]=uigetfile({'*.tif'},'Choose unirradiated film scan');
if (~strcmp(class(ifilename1),'double'))
    flag = 0;
    [I1,delta] = fct_read_tif_image(fct_makecleanfilename(ipathname1,ifilename1),'All');    
    figure('NumberTitle','off','Name','Background');
    h = fct_display(I1,delta);
    ans = questdlg('Do you want to zoom to a smaller ROI?','Image','Yes','No','Yes') ;
    if strcmp(ans,'Yes')
        figure(h);
        I1 = imcrop;
        h = fct_display(I1,delta);
    end
    figure(h);
    [nlines,ncols,ncolors] = size(I1);
    [xgrid1,ygrid1] = fct_gridindextopos(nlines,ncols,delta);
    X = []; Y = [];
    for i=1:4
        [z_zoom,rectoff] = imcrop;
        x1 = delta*round(rectoff(1)/delta);
        x2 = delta*round((rectoff(1)+rectoff(3))/delta);
        y1 = delta*round(rectoff(2)/delta);
        y2 = delta*round((rectoff(2)+rectoff(4))/delta);
        offset = [x1 + x2  y1 + y2]/2;
        figure('NumberTitle','off','Name','Zoom');
        hroi = fct_display(z_zoom,delta);
        clear x y;
        [x,y] = ginput;
        close(hroi);
        if length(x)>1
            x = x(end); y = y(end);
        end
        x = x + offset(1); y = y + offset(2);        
        j = find(abs(xgrid1-x)==min(abs(xgrid1-x)));
        i = find(abs(ygrid1-y)==min(abs(ygrid1-y)));
        X = cat(1,X,j); Y = cat(1,Y,i);
        %check if I aim correctly
%         I1((j-5):(j+5),(i-5):(i+5),:) = 0;
%         h = fct_display(I1,delta);
    end    
    close(h);
    i_bck = Y;
    j_bck = X;
else
    flag = -1;
end

if flag==0

[ifilename1,ipathname1]=uigetfile({'*.tif'},'Choose unirradiated film scan');
    if (~strcmp(class(ifilename1),'double'))
        [I2,delta2] = fct_read_tif_image(fct_makecleanfilename(ipathname1,ifilename1),'All');    
        if delta2~=delta
            errordlg('Image resolutions do not match.');
            flag = 1;
        else
            figure('NumberTitle','off','Name','Irradiated film');
            h = fct_display(I2,delta);
            ans = questdlg('Do you want to zoom to a smaller ROI?','Image','Yes','No','Yes') ;
            if strcmp(ans,'Yes')
                figure(h);
                I2 = imcrop;
                h = fct_display(I2,delta);
            end
            figure(h);
            [nlines,ncols,ncolors] = size(I2);
            [xgrid2,ygrid2] = fct_gridindextopos(nlines,ncols,delta);
            X = []; Y = [];
            for i=1:4
                [z_zoom,rectoff] = imcrop;
                x1 = delta*round(rectoff(1)/delta);
                x2 = delta*round((rectoff(1)+rectoff(3))/delta);
                y1 = delta*round(rectoff(2)/delta);
                y2 = delta*round((rectoff(2)+rectoff(4))/delta);
                offset = [x1 + x2  y1 + y2]/2;
                figure('NumberTitle','off','Name','Zoom');
                hroi = fct_display(z_zoom,delta);
                clear x y;
                [x,y] = ginput;
                close(hroi);
                if length(x)>1
                    x = x(end); y = y(end);
                end
                x = x + offset(1); y = y + offset(2);        
                j = find(abs(xgrid2-x)==min(abs(xgrid2-x)));
                i = find(abs(ygrid2-y)==min(abs(ygrid2-y)));
                X = cat(1,X,j); Y = cat(1,Y,i);
                %check if I aim correctly
        %         I2((j-5):(j+5),(i-5):(i+5),:) = 0;
        %         h = fct_display(I2,delta);
            end
            close(h);
            i_irr = Y;
            j_irr = X;
        end
    else
        flag = 1;
    end
end

if flag==0

    %simplistic registration approach without rotation  
    itmp = sort(i_bck); jtmp = sort(j_bck);
    i1 = floor((itmp(1)+itmp(2))/2):ceil((itmp(3)+itmp(4))/2);
    j1 = floor((jtmp(1)+jtmp(2))/2):ceil((jtmp(3)+jtmp(4))/2);
    i1 = min(i_bck):max(i_bck); j1 = min(j_bck):max(j_bck);    
    Ibck = I1(i1,j1,1);
    Ibck = cat(3,Ibck,I1(i1,j1,2));
    Ibck = cat(3,Ibck,I1(i1,j1,3));
%     Ibck = -log10(max(double(I1(i1,j1,1)),1)./(2^(16)-1));
%     Ibck = cat(3,Ibck,-log10(max(double(I1(i1,j1,2)),1)./(2^(16)-1)));
%     Ibck = cat(3,Ibck,-log10(max(double(I1(i1,j1,3)),1)./(2^(16)-1)));

    itmp = sort(i_irr); jtmp = sort(j_irr);
    i2 = floor((itmp(1)+itmp(2))/2):ceil((itmp(3)+itmp(4))/2);
    j2 = floor((jtmp(1)+jtmp(2))/2):ceil((jtmp(3)+jtmp(4))/2);
    i2 = min(i_irr):max(i_irr); j2 = min(j_irr):max(j_irr);
    Iirr = I2(i2,j2,1);
    Iirr = cat(3,Iirr,I2(i2,j2,2));
    Iirr = cat(3,Iirr,I2(i2,j2,3));
%     Iirr = -log10(max(double(I2(i2,j2,1)),1)./(2^(16)-1));
%     Iirr = cat(3,Iirr,-log10(max(double(I2(i2,j2,2)),1)./(2^(16)-1)));
%     Iirr = cat(3,Iirr,-log10(max(double(I2(i2,j2,3)),1)./(2^(16)-1)));

    [l1,c1,~] = size(Ibck);
    [l2,c2,~] = size(Iirr);

    l = min(l1,l2); c = min(c1,c2);
    off1 = floor((l1-l)/2); l1 = off1+(1:l);
    off2 = floor((l2-l)/2); l2 = off2+(1:l);
    off3 = floor((c1-c)/2); c1 = off3+(1:c);
    off4 = floor((c2-c)/2); c2 = off4+(1:c);    
    
    Ibck = Ibck(l1,c1,:);
    Iirr = Iirr(l2,c2,:);

    h1 = figure;
    h1 = fct_display(Ibck,delta);
    h2 = figure;
    h2 = fct_display(Iirr,delta);



    tmp = date;
    tmp = [tmp '_film_'];
    [ifilename,ipathname]=uiputfile({'*.mat'},'Regestered net image to save',tmp);
    ACTUALDIR = cd;
    if ifilename==0
    else
        filename = fct_makecleanfilename(ipathname,ifilename);        
        save(filename,'Ibck','Iirr','delta');
    end

    err = 1;
elseif flag==1    
    errordlg('Something went wrong. One of the images were not read or the resolutions do not match.')
end