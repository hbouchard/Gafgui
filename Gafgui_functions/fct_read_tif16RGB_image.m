% --------------------------------------------------------------------
function [IMAGE,Resolution,BITS,channel] = fct_read_tif16RGB_image(Filename,signal)

%  fct_read_tif16RGB_image   Read RIT .tif images.
%  The variable signal can be 'Ask', 'Red', Green', 'Blue', 'All'
% IMAGE: the OD value of the film corresponding to the signal
% for exemple, if signal is 'Red/Blue' the value is OD(red) divided by OD(blue)
% Pixel values are converted to double
% Resolution: resolution of the image in cm
% BITS: number of bits for each pixel (ex. 8, 16, 32, ...)
% channel: return the code value for the channel used (1 for red, 2- green, 3- blue,
% 4- white(average of the 3 channels, 5- all channels are kept). For more
% details see function fct_colortochannel
info=imfinfo(Filename)
FORMAT=info.Format
BITS=info.BitsPerSample(1)
NC=size(info.BitsPerSample,2)
DPIX=info.XResolution
DPIY=info.YResolution
COLORTYPE=info.PhotometricInterpretation

if (strcmp(FORMAT,'tif')&&(BITS==16)&&(NC==3)&&(strcmp(COLORTYPE,'RGB'))&&(DPIX==DPIY))
    Resolution=2.54/DPIX;
    IMAGE = imread(Filename);
else
    IMAGE = 0;
    Resolution = 0;
    BITS = 0;
    channel = 0;
end
if NC==3
    if strcmp(signal,'Ask')
        button = questdlg('Which signal would you like to analyze?','Signal','Red','Green','Blue','Red');
    else
        button = signal;
    end
    if strcmp(button,'Red')
        IMAGE=double(IMAGE(:,:,1));
        channel=1;
    elseif strcmp(button,'Green')
        IMAGE=double(IMAGE(:,:,2));
        channel=2;
    elseif strcmp(button,'Blue')
        IMAGE=double(IMAGE(:,:,3));
        channel=3;
    elseif strcmp(button,'All') 
        IMAGE = double(IMAGE);
        channel=5;
    end
else
    IMAGE=double(IMAGE(:,:,1));
    channel=4;
end