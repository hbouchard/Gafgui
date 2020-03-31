% function err = fct_CalibrateSignal(handles);

clear all;
clear functions;
close all;
clc;
fct_AddGafguiFctPath();
handles = fct_initGafgui();

err = 1;

[ifilename,ipathname]=uigetfile({'*.tif'},'Film to calibrate');

if ~strcmp(class(ifilename),'double')
    
    [Im,DPIX,BITS,channel] = fct_read_tif_image(fct_makecleanfilename(ipathname,ifilename),'All');
    Im = double(Im)/(2^BITS-1);
    
    [cfilename,cpathname]=uigetfile({'*.od'},'OD calibration file');    
    if (~strcmp(class(ifilename),'double'))
        Filename = fct_makecleanfilename(cpathname,cfilename);
        file = fopen(Filename,'r');
        [I,Sr,Sg,Sb] = fct_ReadODfile(file);
        fclose(file);
        flag = 2;
    else
        errordlg('Wrong OD calibration file. Signal will not be corrected.')
        I = [0 1];
        Sr = [0 1];
        Sg = [0 1];
        Sb = [0 1];
    end
    Im(:,:,1) = uint16((2^BITS-1)*interp1(Sr,I,Im(:,:,1),'pchip'));
    Im(:,:,2) = uint16((2^BITS-1)*interp1(Sg,I,Im(:,:,2),'pchip'));
    Im(:,:,3) = uint16((2^BITS-1)*interp1(Sb,I,Im(:,:,3),'pchip'));
    
    tfilename = [ifilename(1:(end-4)) '_cal.tif'];
    
    [ofilename,opathname] = uiputfile({'*.tif'},'Save average film',string(tfilename));
    if ofilename==0
    else
        ofilename = fct_makecleanfilename(opathname,ofilename);
        imwrite(Im,ofilename,'tif','Resolution',round(2.54/DPIX));
    end

end
