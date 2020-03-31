% --------------------------------------------------------------------
function [IMAGE,Resolution,BITS,channel] = fct_read_tif_image(Filename,signal)

%  fct_read_tif16RGB_image   Read .tif images.
%  The variable signal can be 'Red', Green', 'Blue', 'All'
% Id it isn't one of these, it will ask
% IMAGE: the OD value of the film corresponding to the signal
% for exemple, if signal is 'Red/Blue' the value is OD(red) divided by OD(blue)
% Pixel values are converted to double
% Resolution: resolution of the image in cm
% BITS: number of bits for each pixel (ex. 8, 16, 32, ...)
% channel: return the code value for the channel used (1 for red, 2- green, 3- blue,
% 4- white(average of the 3 channels, 5- all channels are kept). For more
% details see function fct_colortochannel
info = imfinfo(Filename)
FORMAT = info.Format
BITS = info.BitsPerSample(1)
NC = size(info.BitsPerSample,2)
DPIX = info.XResolution
DPIY = info.YResolution
COLORTYPE = info.PhotometricInterpretation

if (strcmp(FORMAT,'tif')&&(NC==3)&&(strcmp(COLORTYPE,'RGB'))&&(DPIX==DPIY))
    Resolution=2.54/DPIX;
    IMAGE = imread(Filename);
elseif (strcmp(FORMAT,'tif')&&(NC==1)&&(strcmp(COLORTYPE,'BlackIsZero'))&&(DPIX==DPIY))
    Resolution=2.54/DPIX;
    IMAGE = imread(Filename);
else
    error('Wrong file format.');
end
if NC==3
    if (~strcmp(signal,'Red'))&&(~strcmp(signal,'Green'))&&(~strcmp(signal,'Blue'))&&(~strcmp(signal,'All'))
        button = questdlg('Which signal would you like to analyze?','Signal','Red','Green','Blue','Red');
        if length(button) == 0
            button = 'Red';
        end
%         str{1} = 'Red';
%         str{2} = 'Green';
%         str{3} = 'Blue';
%         str{4} = 'Multi';
%         [button,ok] = listdlg('Name','Channel','ListString',str);
%         if (ok==1)&&(length(button)==1)
%             button = str{button};
%         else button=[];
%         end
    else
        button = signal;
    end
    if strcmp(button,'Red')
        IMAGE = IMAGE(:,:,1);
        channel = 1;
    elseif strcmp(button,'Green')
        IMAGE = IMAGE(:,:,2);
        channel = 2;
    elseif strcmp(button,'Blue')
        IMAGE = IMAGE(:,:,3);
        channel = 3;
    elseif strcmp(button,'All') 
        IMAGE = IMAGE;
        channel=5;
    end
elseif NC==1
    IMAGE = IMAGE;
    channel = 4;
else
    error('Number of channels is neither 1 or 3.');
end