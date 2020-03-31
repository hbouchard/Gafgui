function handles = fct_MultiCorrect(handles);
% % 
% clc;
% clear all;
% close all;
% fct_AddGafguiFctPath();
% handles = fct_initGafgui();
% close;

%%%%%%%%%%%%%%%%%%%%%%%%%
%Open an image to correct
%%%%%%%%%%%%%%%%%%%%%%%%%

[ifilename,ipathname]=uigetfile({'*.tif'},'Image to correct');
if (~strcmp(class(ifilename),'double'))
    flag = 2;
else
    flag = 0;
end
if flag == 2
    ACTUALDIR=cd;
    % details see function fct_colortochannel
    Filename = fct_makecleanfilename(ipathname,ifilename)
    info=imfinfo(Filename);
    FORMAT=info.Format;
    BITS=info.BitsPerSample(1);
    NC=size(info.BitsPerSample,2);
    DPIX=info.XResolution;
    DPIY=info.YResolution;
    COLORTYPE=info.PhotometricInterpretation;

    if (strcmp(FORMAT,'tif')&&(BITS==16)&&(NC==3)&&(strcmp(COLORTYPE,'RGB'))&&(DPIX==DPIY))
        delta = 2.54/DPIX;
        IMAGE = imread(Filename);
%         IMAGE(:,:,1) = (IMAGE(:,:,1)); IMAGE(:,:,2) = (IMAGE(:,:,2));  IMAGE(:,:,3) = (IMAGE(:,:,3));
    else
        error('Wrong file format.');
    end

    [ofilename,opathname]=uigetfile({'*.3ch'},['Load mutlichannel correction method']);
    if strcmp(class(ofilename),'double')
    else
        fname = fct_makecleanfilename(opathname,ofilename);
        %new method
        
        [Rot,rrange,grange,brange,trange,err]  = fct_ReadMultiCorrection(fname);
        if err
            errordlg('File not read. Something went wrong.');
        else
            %%%%%%%%%%%%%%%%%%
            %prepare filter to reduce noize
            clear str;

            %button = questdlg('Do you want to filter the thickness correction?','Filter','Yes','No','Yes') ;
            button = 'No';
            if strcmp(button,'Yes')
                filterimage = 1;
            else
                filterimage = 0;
            end
            if filterimage
                kconv = max(1,floor(0.1/delta));
                kconv = 2*floor(kconv/2)+1;
            else
                kconv = 1;
            end
            filter = ones(kconv,kconv)/sum(sum(ones(kconv,kconv)));
            dummy = filter*0;
            dummy(floor(kconv/2)+1,floor(kconv/2)+1) = 1;

            %This creates a 2D flag being 1 if valid range and zero
            %otherwise.
            R = double(IMAGE(:,:,1));
            G = double(IMAGE(:,:,2));
            B = double(IMAGE(:,:,3)); 
            FLAG = ones(size(R,1),size(R,2));
            FLAG = FLAG.*double(logical(R>=rrange(1))).*double(logical(R<=rrange(2)));
            FLAG = FLAG.*double(logical(G>=grange(1))).*double(logical(G<=grange(2)));
            FLAG = FLAG.*double(logical(B>=brange(1))).*double(logical(B<=brange(2)));
            
            %now the OD
            R = -log10(double(IMAGE(:,:,1))/65535); 
            G = -log10(double(IMAGE(:,:,2))/65535);  
            B = -log10(double(IMAGE(:,:,3))/65535);

            X = Rot(1,1)*R + Rot(1,2)*G + Rot(1,3)*B; 
            Z = Rot(3,1)*R + Rot(3,2)*G + Rot(3,3)*B; 

            X = conv2(X,dummy,'valid');
            Z = conv2(Z,filter,'valid');
            Z = Z/mean(Z(:));
            %we want to avoid Z=0
            [lin,col] = size(Z);
            Z = reshape(Z,lin*col,1);
            Z(find(Z==0)) = 1;
            Z = reshape(Z,lin,col);
            THETA = X./Z;
            FLAG = FLAG.*double(logical(THETA>=trange(1))).*double(logical(THETA<=trange(2)));
            IMG = uint16(FLAG.*(65535*10.^(-THETA)) + (1-FLAG).*(65535*10.^(-R)));
            %
            if 0%to show the difference between corrected and uncorrected
                defaultname = sprintf('%s_raw.tif',ifilename(1:(length(ifilename)-4)));
                [ofilename,opathname] = uiputfile({'*.tif'},'.tif image to export without correction',defaultname);
                if strcmp(class(ofilename),'double')
                else
                    filename = fct_makecleanfilename(opathname,ofilename);
                    Imraw = uint16(65535*10.^(-R));
    %                 Imraw(:,:,1) = uint16(65535*10.^(-R));
    %                 Imraw(:,:,2) = uint16(65535*10.^(-G));
    %                 Imraw(:,:,3) = uint16(65535*10.^(-B));
                    imwrite(Imraw,filename,'tif','Resolution',floor(2.54/delta));
                end
            end
            defaultname = sprintf('%s_multi.tif',ifilename(1:(length(ifilename)-4)));
            [ofilename,opathname] = uiputfile({'*.tif'},'.tif image to export with correction',defaultname);
            if strcmp(class(ofilename),'double')
            else
                filename = fct_makecleanfilename(opathname,ofilename);
                Imcorr = IMG;
%                 Imcorr(:,:,1) = uint16(65535*10.^(-THETA));
%                 Imcorr(:,:,2) = uint16(65535*10.^(-THETA));
%                 Imcorr(:,:,3) = uint16(65535*10.^(-THETA));
                imwrite(Imcorr,filename,'tif','Resolution',floor(2.54/delta));
            end
        end

    end
end
