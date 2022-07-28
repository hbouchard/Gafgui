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
        
%         [Rot,rrange,grange,brange,trange,znorm,err]  = fct_ReadMultiCorrection(fname);
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
                kconv = max(1,floor(0.1/delta));%1 cm
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
            %HB: 11 July 2022: bug fix to assure we can use film for which OD values are too close to range ends 
            %FLAG = FLAG.*double(logical(R>=rrange(1))).*double(logical(R<=rrange(2)));
            %FLAG = FLAG.*double(logical(G>=grange(1))).*double(logical(G<=grange(2)));
            %FLAG = FLAG.*double(logical(B>=brange(1))).*double(logical(B<=brange(2)));
            %HB: the trick is to reduce the distance of limits from 0 and 65535 by half
            SMIN = min([rrange(:)' grange(:)' brange(:)']);
            SMAX = max([rrange(:)' grange(:)' brange(:)']);
            SMIN = floor(SMIN - SMIN/2);
            SMAX = ceil(SMAX + (65635-SMAX)/2);
            FLAG = FLAG.*double(logical(R>=SMIN)).*double(logical(R<=SMAX));
            FLAG = FLAG.*double(logical(G>=SMIN)).*double(logical(G<=SMAX));
            FLAG = FLAG.*double(logical(B>=SMIN)).*double(logical(B<=SMAX));
            %now the OD
            R = -log10(double(IMAGE(:,:,1))/65535); 
            G = -log10(double(IMAGE(:,:,2))/65535);  
            B = -log10(double(IMAGE(:,:,3))/65535);
            %the eigencolors
            X = Rot(1,1)*R + Rot(1,2)*G + Rot(1,3)*B; 
            Z = Rot(3,1)*R + Rot(3,2)*G + Rot(3,3)*B; 
            %convolution of EC3 with a small filter whilst EC1 remain exact
            X = conv2(X,dummy,'valid');
            Z = conv2(Z,filter,'valid');
            FLAG = conv2(FLAG,dummy,'valid');
            S = conv2(R,dummy,'valid');
            %Bug William 11 juillet 2022??
            %Z = Z/mean(Z(:));%HB says: OH MY GOD NO!
            %HB 11 July 2022: this is crucial because we need to keep theta
            %within the range as we save the image into tiff format, hence use
            %I = 65535*10.^(-THETA). One alternative would be to ask the
            %user to select a piece of film (exempt of marks or tape) and take the average EC3 for
            %normalization. Until then, we pass the normalization information from the characterization of the correction 
            %Z = Z/znorm;
            %HB 12 July: we will substract by the minimum value of trange
            %instead and add 0.1
            %we want to avoid Z=0
            [lin,col] = size(Z);
            Z = reshape(Z,lin*col,1);
            Z(find(Z==0)) = 1;
            Z = reshape(Z,lin,col);
            THETA = X./Z;
            %Bug William 11 juillet 2022
            %FLAG = FLAG.*double(logical(THETA>=trange(1))).*double(logical(THETA<=trange(2)));
            %HB 12 July: this is the max OD for which there will be a numerical error
            %of 0.0001 on OD via conversion into unint16
            maxOD = log10(10^0.0001-1)+log10(65535);
            minOD = 0.1;
            %HB 12 July this is the new way to assure THETA can be
            %converted into a 16 bit tiff image without loosing information
            %hence we map THETA from trange to [minOD,maxOD]
            THETA = (THETA - min(trange))/(max(trange)- min(trange))*(maxOD-minOD) + minOD;
            THETA = max(min(THETA,log10(65535)),0);
            %THETA = max(min(THETA,maxOD*2),0);
            % Now we use the valid part only for which RGB falls into range
            % which assures there is actually film (or something else) on the scanner bed
            IMG = uint16(FLAG.*(65535*10.^(-THETA)) + (1-FLAG).*(65535*10.^(-S)));
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
