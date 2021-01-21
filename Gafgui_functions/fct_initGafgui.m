function handles = fct_initGafgui(varargin)

%Initialization of global variables (defined only once)
handles = [];
if numel(varargin)==0   
    handles.defaultcolor = 'Multi';
elseif numel(varargin)==1
    handles.defaultcolor = varargin{1};    
end

%Initialization of main variables
handles.CCDres = 0.0010583; %in cm
handles.dataini = 0; %data type initialized  1=signal, 2=OD, 3=net OD, 4=DOSE
handles.SATMINini = 0; %minimum saturation value initialized
handles.SATMAXini = 0; %maximum saturation value initialized
handles.BCKGRNDini = []; % OD value of the unexposed film initialized
handles.NpixBCKGRNDini = 0; %number of pixels in the background value initialized (if we use large ROIs, this is no longer necessary since Bouchard et al 2021 include the incertainty in sig0)
handles.DOSEREFini = []; %reference dose value for normalization initialized
handles.NpixDOSEREFini = []; %number of pixels in the reference dose value for normalization initialized
handles.DELTAini = 0; %resolution initialized
handles.channelini = fct_colortochannel(handles.defaultcolor); %the channel initialized 1=red, 2=green, 3=blue
handles.filterini = 0; %filter size being applied initialized
handles.ORIGINini = [0 0]'; %origin position initialized
handles.Z = []; %the image: Z is the image as loaded/imported
handles.cfilenameini = '';% the calibration file name initialized
handles.BITSini = 0; % bits number of the image initialized ex. 8, 16, 32, ...
handles.data = handles.dataini;%the initial environement variables
handles.SATMIN = handles.SATMINini;%the minimum saturation value
handles.SATMAX = handles.SATMAXini;%the maximum saturation value
handles.BCKGRND = handles.BCKGRNDini;% OD value of the unexposed film
handles.NpixBCKGRND = handles.NpixBCKGRNDini;%number of pixels in the background value (if we use large ROIs, this is no longer necessary since Bouchard et al 2021 include the incertainty in sig0)
handles.DOSEREF = handles.DOSEREFini;%reference dose value for normalization
handles.NpixDOSEREF = handles.NpixDOSEREFini;%number of pixels in the reference dose value for normalization
handles.DELTA = handles.DELTAini; %the resolution
handles.channel = handles.channelini; %the channel 1=red, 2=green, 3=blue
handles.filter = handles.filterini; %filter size being applied 
handles.ORIGIN = handles.ORIGINini;%the position of origin
handles.z = handles.Z;  %the image: z is the image as it is currently
handles.cfilename = handles.cfilenameini;% the calibration file name
%Initialization of other variables
handles.color='jet';% the colormap
handles.gridflag = 0;% grin on/off flag
handles.gridres = 1;%the grid resolution in cm
handles.ifilename = ''; %the input filename as imported
handles.H = gcf; %current handles
handles.menus = guihandles(handles.H); %the menus
handles = fct_updatedisplay(handles); %the whole se of variables
