% function err = fct_SubtractBackground(handles);
% 
clear all;
clear functions;
close all;
clc;
fct_AddGafguiFctPath();
handles = fct_initGafgui();

flag = 0;
[ifilename1,ipathname1]=uigetfile({'*.tif'},'Choose unirradiated film');
[ifilename2,ipathname2]=uigetfile({'*.tif'},'Choose irradiated film');
if (~strcmp(class(ifilename1),'double'))&&(~strcmp(class(ifilename2),'double'))
    
    flag = 1;
    
    [Im1,RES1,BITS1] = fct_read_tif_image(fct_makecleanfilename(ipathname1,ifilename1),'All');
    [Im2,RES2,BITS2] = fct_read_tif_image(fct_makecleanfilename(ipathname2,ifilename2),'All');
    
    [nlines1,ncols1] = size(Im1(:,:,1));
    [xgrid1,ygrid1] = fct_gridindextopos(nlines1,ncols1,RES1);
    [nlines2,ncols2] = size(Im2(:,:,1));
    [xgrid2,ygrid2] = fct_gridindextopos(nlines2,ncols2,RES2);
    if (RES1~=RES2)||(BITS1~=16)||(BITS2~=16)
       errodlg('Image resoltuion mismatch.');
       flag = 0; 
    end
    if flag==1
        RES = RES1; clear RES1 RES2;
        h = fct_display(Im1,RES);
        h1 = msgbox('Click on all 4 corners of the unirradiated film piece and press enter','Positions');
        uiwait(h1);
        [x1,y1,nx1,ny1] = fct_getpoints(Im1,RES);
        close(h);
        if length(x1)~=4
            errodlg('You need 4 points.');
            flag = 0;
        end
    end
    if flag==1
        h = fct_display(Im2,RES);
        h2 = msgbox('Click on all 4 corners of the irradiated film piece and press enter','Positions');
        uiwait(h2);
        [x2,y2,nx2,ny2] = fct_getpoints(Im2,RES);
        close(h);
        if length(x2)~=4
            errodlg('You need 4 points.');
            flag = 0;
        end
    end
    if flag==1
       jmin1 = fct_postoindex(min(x1),xgrid1); jmax1 = fct_postoindex(max(x1),xgrid1); 
       imin1 = fct_postoindex(min(y1),ygrid1); imax1 = fct_postoindex(max(y1),ygrid1); 
       I1 = Im1(imin1:imax1,jmin1:jmax1,:);
       jmin2 = fct_postoindex(min(x2),xgrid2); jmax2 = fct_postoindex(max(x2),xgrid2); 
       imin2 = fct_postoindex(min(y2),ygrid2); imax2 = fct_postoindex(max(y2),ygrid2); 
       I2 = Im2(imin2:imax2,jmin2:jmax2,:);
       [m1,n1,ch] = size(I1);
       [m2,n2,ch] = size(I2);
       m = min(m1,m2); n = min(n1,n2);
       I1 = I1(floor((m1-m)/2)+(1:m),floor((n1-n)/2)+(1:n),:);
       I2 = I2(floor((m2-m)/2)+(1:m),floor((n2-n)/2)+(1:n),:);
       clear I;
       netR = log10(65535./double(I2(:,:,1)))-log10(65535./double(I1(:,:,1)));
       netG = log10(65535./double(I2(:,:,2)))-log10(65535./double(I1(:,:,2)));
       netB = log10(65535./double(I2(:,:,3)))-log10(65535./double(I1(:,:,3)));
       I(:,:,1) = uint16(65535*10.^(-netR));
       I(:,:,2) = uint16(65535*10.^(-netG));
       I(:,:,3) = uint16(65535*10.^(-netB));
       figure;
       subplot(1,3,1);
       imagesc(I1);
       set(gca,'DataAspectRatio',[1 1 1]);
       subplot(1,3,2);
       imagesc(I2);
       set(gca,'DataAspectRatio',[1 1 1]);
       subplot(1,3,3);
       imagesc(I);
       set(gca,'DataAspectRatio',[1 1 1]);
       impixelinfo;
    end
    if flag==1
        ifilename2 = ifilename2(:)';
        tfilename = [ifilename2(1:(end-4)) '_net.tif']
        [ofilename,opathname] = uiputfile({'*.tif'},'Save average film',string(tfilename));
        if ofilename==0
        else
            ofilename = fct_makecleanfilename(opathname,ofilename);
            imwrite(I,ofilename,'tif','Resolution',round(2.54/RES));
        end
    end
end
   