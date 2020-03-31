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
handles.dataini = 0;
handles.SATMINini = 0;
handles.SATMAXini = 0;
handles.BCKGRNDini = []; % OD value of the unexposed film
handles.NpixBCKGRNDini = 0;
handles.DOSEREFini = [];
handles.NpixDOSEREFini = [];
handles.DELTAini = 0;
handles.channelini = fct_colortochannel(handles.defaultcolor);
handles.filterini = 0;
handles.ORIGINini = [0 0]';
handles.Z = [];
handles.cfilenameini = '';
handles.BITSini = 0; % bits number of the image ex. 8, 16, 32, ...
handles.data = handles.dataini;
handles.SATMIN = handles.SATMINini;
handles.SATMAX = handles.SATMAXini;
handles.BCKGRND = handles.BCKGRNDini;
handles.NpixBCKGRND = handles.NpixBCKGRNDini;
handles.DOSEREF = handles.DOSEREFini;
handles.NpixDOSEREF = handles.NpixDOSEREFini;
handles.DELTA = handles.DELTAini;
handles.channel = handles.channelini;
handles.filter = handles.filterini;
handles.ORIGIN = handles.ORIGINini;
handles.z = handles.Z;
handles.cfilename = handles.cfilenameini;
%Initialization of other variables
handles.color='jet';
handles.gridflag = 0;
handles.gridres = 1;
handles.ifilename = '';
handles.H = gcf;
handles.menus = guihandles(handles.H);
handles = fct_updatedisplay(handles);
