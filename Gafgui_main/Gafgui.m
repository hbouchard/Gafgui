%     Gafgui version 4.0
%     Copyright (C) 2007-2020  by Hugo Bouchard. 
%     My contact is h.bouchard@umontreal.ca
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%     If you need information on the science behind this code, or if you
%     publish results obtained with that program, please see or cite this
%     reference:
%
%     H. Bouchard, F. Lacroix, G. Beaudoin, J. Carrier, I. Kawrakow, On the
%     characterization and uncertainty analysis of radiochromic film
%     dosimetry, Med. Phys. Volume 36, Issue 6, pp. 1931-1946 (June 2009)
%

function varargout = Gafgui(varargin)
% Gafgui M-file for Gafgui.fig
%      Gafgui, by itself, creates a new Gafgui or raises the existing
%      singleton*.
%
%      H = Gafgui returns the handle to a new Gafgui or the handle to
%      the existing singleton*.
%
%      Gafgui('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in Gafgui.M with the given input arguments.
%
%      Gafgui('Property','Value',...) creates a new Gafgui or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Gafgui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Gafgui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Gafgui

% Last Modified by GUIDE v2.5 11-Feb-2020 18:50:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Gafgui_OpeningFcn, ...
    'gui_OutputFcn',  @Gafgui_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
%addpath('H:\dsp\Radio-oncologie\commun\Physique radio-onco\MATLAB\SPRO_functions');
%addpath('C:\Program Files\MATLAB\R2007a\toolbox\spm5')
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Gafgui is made visible.
function Gafgui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Gafgui (see VARARGIN)

format long;
warning off;
clc;


%%%%%%%%%%%%%
%This add the functions path
cpu = computer;
if strcmp(cpu,'PCWIN')||strcmp(cpu,'PCWIN64')
    c = '\'; % dos
else
    c = '/'; % mac/linux/unix
    %In the case someone uses a SUN, I don't know if this is correct
end
pth = cd;
[lin,col] = size(pth);
if lin>col
    pth = pth';
end
l = length(pth);
if pth(end)==c
    pth = pth(1:end-1);
end
k = max(find(pth==c))-1;
pth = pth(1:k);
pth = [pth c 'Gafgui_functions'];
addpath(pth);
%%%%%%%%%%%%%%

%Agree to conditions
if fct_version()==1
    close all;
    rmpath(pth);
else
    handles = fct_initGafgui();
    % Choose default command line output for Gafgui
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);

    % UIWAIT makes Gafgui wait for user response (see UIRESUME)
    % uiwait(handles.figure1);

    % --- Outputs from this function are returned to the command line.
end

function varargout = Gafgui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isfield(handles,'output');
    varargout{1} = handles.output;
else
    varargout{1} = [];%I have to put this in case the gui doesn't pass the licence agreement
end

% --------------------------------------------------------------------
function MAINFIGURE_CreateFcn(hObject, eventdata, handles, varargin)

% --------------------------------------------------------------------
function MAINFIGURE_ButtonDownFcn(hObject, eventdata, handles)

% --------------------------------------------------------------------
function MENUFILE_Callback(hObject, eventdata, handles)
% hObject    handle to MENUFILE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MENUIMPORTIMAGE_Callback(hObject, eventdata, handles)
% hObject    handle to MENUIMPORTIMAGE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUTOOLS_Callback(hObject, eventdata, handles)
% hObject    handle to MENUTOOLS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MEAN_MEANDOSE_ORIGIN_Callback(hObject, eventdata, handles)
% hObject    handle to MEAN_MEANDOSE_ORIGIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUROT_Callback(hObject, eventdata, handles)
% hObject    handle to MENUROT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUFILTERIM_Callback(hObject, eventdata, handles)
% hObject    handle to FILTERIM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUROI_Callback(hObject, eventdata, handles)
% hObject    handle to MENUROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUANALYZE_Callback(hObject, eventdata, handles)
% hObject    handle to MENUANALYZE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUDIST2D_Callback(hObject, eventdata, handles)
% hObject    handle to ANALYDIST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENU_MEANDOSE_Callback(hObject, eventdata, handles)
% hObject    handle to MENU_MEANDOSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUDONNE_Callback(hObject, eventdata, handles)
% hObject    handle to MENUDONNE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUAFF_Callback(hObject, eventdata, handles)
% hObject    handle to MENUAFF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUABOUT_Callback(hObject, eventdata, handles)
% hObject    handle to MENUABOUT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUCALIBRATION_Callback(hObject, eventdata, handles)
% hObject    handle to MENUCALIBRATION (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MENURP_Callback(hObject, eventdata, handles)
% hObject    handle to MENURP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------

function MENUPROFIL_Callback(hObject, eventdata, handles)
% hObject    handle to MENUPROFIL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUPREP_Callback(hObject, eventdata, handles)
% hObject    handle to MENUPREP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENU_SATURATE_Callback(hObject, eventdata, handles)
% hObject    handle to MENU_SATURATE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENU_SAVE_IMAGE_Callback(hObject, eventdata, handles)
% hObject    handle to MENU_SAVE_IMAGE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUREGISTRATION_Callback(hObject, eventdata, handles)
% hObject    handle to MENUREGISTRATION (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUPROFILE_Callback(hObject, eventdata, handles)
% hObject    handle to MENUPROFILE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENU_RP_Callback(hObject, eventdata, handles)
% hObject    handle to MENURP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENU2DDIST_Callback(hObject, eventdata, handles)
% hObject    handle to MENU2DDIST (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUEXPORT_Callback(hObject, eventdata, handles)
% hObject    handle to MENUEXPORT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUHOMOG_Callback(hObject, eventdata, handles)
% hObject    handle to MENUHOMOG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUCHARACT_Callback(hObject, eventdata, handles)
% hObject    handle to MENUCHARACT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENUPARAM_Callback(hObject, eventdata, handles)
% hObject    handle to MENUPARAM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MENU_VARANALYSIS_Callback(hObject, eventdata, handles)
% hObject    handle to MENU_VARANALYSIS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function MULTICHANNEL_CHARACTERIZE_Callback(hObject, eventdata, handles)
% hObject    handle to MULTICHANNEL_CHARACTERIZE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%err = fct_MultichCharactLinearNew(handles);
handles = fct_CreateMultichannelCorrection(handles)
% --------------------------------------------------------------------
function MULTICHANNEL_CORRECT_Callback(hObject, eventdata, handles)
% hObject    handle to MULTICHANNEL_CORRECT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles = fct_MultichCorrectLinear(handles);
% handles = fct_MultichCorrectRatioGreenBlue(handles);
handles = fct_MultiCorrect(handles);

% --------------------------------------------------------------------
function QUIT_Callback(hObject, eventdata, handles)
% hObject    handle to QUIT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clear handles;
clear all;
close all;

% --------------------------------------------------------------------
function RESTOREIM_Callback(hObject, eventdata, handles)
% hObject    handle to RESTOREIM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    handles.data = handles.dataini;
    handles.SATMIN = handles.SATMINini;
    handles.SATMAX = handles.SATMAXini;
    %handles.BCKGRND = handles.BCKGRNDini;
    %handles.NpixBCKGRND = handles.NpixBCKGRNDini;
    %handles.DOSEREF = handles.DOSEREFini;
    %handles.NpixDOSEREF = handles.NpixDOSEREFini;
    handles.DELTA = handles.DELTAini;
    handles.channel = handles.channelini;
    handles.filter = handles.filterini;
    handles.ORIGIN = handles.ORIGINini;
    handles.z = handles.Z;
    handles.cfilename = handles.cfilenameini;
    %Other variables
    if handles.channel==4
        handles.color = 'gray';
    else
        handles.color ='jet';
    end
    handles.gridflag = 0;
    handles.gridres = fct_autoruler(handles.z,handles.DELTA);
    %handles.ifilename = fct_makecleanfilename(ipathname,ifilename);
    handles = fct_updatedisplay(handles);;
    % Update handles structure
    guidata(hObject, handles);
    
end

% --------------------------------------------------------------------
function ROTRIGHT_Callback(hObject, eventdata, handles)
% hObject    handle to ROTRIGHT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS
    z = fliplr(handles.z');
    %CREATION OF NEW VARIABLES
    handles.z = z;
    %UPDATE DISPLAY
    handles = fct_updatedisplay(handles);;
    % Update handles structure
    guidata(hObject, handles);
    
end

% --------------------------------------------------------------------
function ROTLEFT_Callback(hObject, eventdata, handles)
% hObject    handle to ROTLEFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS
    z=flipud(handles.z');
    %CREATION OF NEW VARIABLES
    handles.z=z;
    %UPDATE DISPLAY
    handles = fct_updatedisplay(handles);;
    % Update handles structure
    guidata(hObject, handles);
    
end


% --------------------------------------------------------------------
function FOLD_IMAGE_Callback(hObject, eventdata, handles)
% hObject    handle to FOLD_IMAGE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if fct_isthereanimage(handles)
    
    str{1} = 'Vertical';
    str{2} = 'Horizontal'
    str{3} = '360';
    [button,ok] = listdlg('Name','Direction','ListString',str);
    
    if (ok==1)&&(length(button)==1)
        button = str{button};
        if(strcmp(button,'Vertical'))
            z = (handles.z +  flipud(handles.z))/2;
        elseif(strcmp(button,'Horizontal'))
            z = (handles.z +  fliplr(handles.z))/2;
        elseif(strcmp(button,'360'))
            z = handles.z;
            [M,N] = size(z);
            del = asin(1/(N/2))*(180/3.141592653589793);
            del = max(del,1);
            tmp = sprintf('Rotation: steps of %.1f deg',del);
            hwait = waitbar(0,tmp);
            for angle = 0:del:90
                figure(hwait);
                waitbar(angle/90,hwait);
                z = z + imrotate(handles.z,angle,'bilinear','crop');
            end
            close(hwait)
            z = z/(length(0:del:90));
            z = (z +  fliplr(z))/2;
            z = (z +  flipud(z))/2;
        end
        %CREATION OF NEW VARIABLES
        handles.z=z;
        %UPDATE DISPLAY
        handles = fct_updatedisplay(handles);;
        % Update handles structure
        guidata(hObject, handles);
    end    
end


% --------------------------------------------------------------------
function MIRVERT_Callback(hObject, eventdata, handles)
% hObject    handle to MIRVERT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS
    z= flipud(handles.z);
    %CREATION OF NEW VARIABLES
    handles.z=z;
    %UPDATE DISPLAY
    handles = fct_updatedisplay(handles);;
    % Update handles structure
    guidata(hObject, handles);
    
end

% --------------------------------------------------------------------
function MIRHORIZ_Callback(hObject, eventdata, handles)
% hObject    handle to MIRHORIZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS
    z=fliplr(handles.z);
    %CREATION OF NEW VARIABLES
    handles.z=z;
    %UPDATE DISPLAY
    handles = fct_updatedisplay(handles);;
    % Update handles structure
    guidata(hObject, handles);
    
end

% --------------------------------------------------------------------
function ROTPOINTS_Callback(hObject, eventdata, handles)
% hObject    handle to ROTPOINTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %TAKE POINTS FOR ROTATION
    [x,y,nx,ny] = fct_getpoints(handles.z,handles.DELTA);

    if length(y)<2
    else
        %take last 2 points
        l = length(y);
        x = x(l-1:l);
        y = y(l-1:l);
        %OPERATIONS
        if (y(2)-y(1))==0
        else
            a=(abs(y(2)-y(1))/(y(2)-y(1)))*acos(sqrt((x(2)-x(1))^2)/sqrt((x(2)-x(1))^2+(y(2)-y(1))^2));
            angle = a*180/3.141592653589793;
            z = imrotate(handles.z,angle,'bilinear','crop');
            handles.z=z;
        end
        %UPDATE DISPLAY
        handles = fct_updatedisplay(handles);;
        % Update handles structure
        guidata(hObject, handles);
    end
end

% --------------------------------------------------------------------
function ROIFIXED_Callback(hObject, eventdata, handles)
% hObject    handle to ROIFIXED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS & CREATION OF NEW VARIABLES
    figure(handles.H);
    
    [newz,center,owidth] = fct_getroi(handles.z,handles.DELTA,'Fixed',0);
    handles.z = newz;
    figure(handles.H);
    %UPDATE DISPLAY
    handles.gridres = fct_autoruler(handles.z,handles.DELTA);
    guidata(hObject, handles);
    handles = fct_updatedisplay(handles);;
    % Update handles structure
    guidata(hObject, handles);
    
end

% --------------------------------------------------------------------
function ROISELECT_Callback(hObject, eventdata, handles)
% hObject    handle to ROISELECT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    
    str{1} = 'Free';
    str{2}='Rectangular';
    str{3}='Zoom';
    str{4} = 'Point';
    
    [type,ok] = listdlg('Name','ROI','ListString',str);
    if (ok==1)
        
        type = str{type};
        
        %OPERATIONS & CREATION OF NEW VARIABLES
        figure(handles.H);
        [newz,center,owidth] = fct_getroi(handles.z,handles.DELTA,type,0);
        [m,n] = size(newz);
        if m*n~=0
        handles.z = newz;
            %UPDATE DISPLAY
            res = fct_autoruler(handles.z,handles.DELTA);
            handles.gridres = res;
            guidata(hObject, handles);
            handles = fct_updatedisplay(handles);;
            % Update handles structure
            guidata(hObject, handles);
        end
    end
end

% --------------------------------------------------------------------
function ROIRECT_Callback(hObject, eventdata, handles)
% hObject    handle to ROIRECT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS & CREATION OF NEW VARIABLES
    figure(handles.H);
    [newz,center,owidth] = fct_getroi(handles.z,handles.DELTA,'Rectangular',0);
    handles.z = newz;
    figure(handles.H);
    %UPDATE DISPLAY
    handles.gridres = fct_autoruler(handles.z,handles.DELTA);
    guidata(hObject, handles);
    handles = fct_updatedisplay(handles);;
    % Update handles structure
    guidata(hObject, handles);
    
end

% --------------------------------------------------------------------
% function IMPORT_RIT_Callback(hObject, eventdata, handles)
% hObject    handle to Import_RIT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% [ifilename,ipathname]=uigetfile({'*.rv4'},'Film to measure');
% 
% if ~strcmp(class(ifilename),'double')
%     
%     ACTUALDIR=cd;
%     handles.ifilename = fct_makecleanfilename(ipathname,ifilename);
%     [z,DPIX,sx,sy] = fct_read_rit_image(handles.ifilename);
%     channel = 4;
%     BITS = 16;
%     %Initial main variables
%     handles.dataini = 2;
%     handles.SATMINini = 0;
%     handles.SATMAXini = 0;
%     handles.BCKGRNDini = [];
%     handles.NpixBCKGRNDini = 0;
%     handles.DOSEREFini = [];
%     handles.NpixDOSEREFini = [];
%     handles.DELTAini = DPIX;
%     handles.channelini = channel;
%     handles.filterini = 0;
%     handles.ORIGINini = [0 0 ]';
%     handles.Z = -log10(max(z,1)./(2^BITS-1));
%     handles.cfilenameini = '';
%     handles.BITSini = BITS;
%     %Actual main variables
%     handles.data = handles.dataini;
%     handles.SATMIN = handles.SATMINini;
%     handles.SATMAX = handles.SATMAXini;
%     %handles.BCKGRND = handles.BCKGRNDini;
%     %%handles.NpixBCKGRND = handles.NpixBCKGRNDini;
%     handles.DOSEREF = handles.DOSEREFini;
%     %handles.NpixDOSEREF = handles.NpixDOSEREFini;
%     handles.DELTA = handles.DELTAini;
%     handles.channel = handles.channelini;
%     handles.filter = handles.filterini;
%     handles.ORIGIN = handles.ORIGINini;
%     handles.z = handles.Z;
%     handles.cfilename = handles.cfilenameini;
%     %Other variables
%     handles.color = 'gray';
%     handles.gridflag = 0;
%     handles.gridres = fct_autoruler(handles.z,handles.DELTA);
%     handles.ifilename = fct_makecleanfilename(ipathname,ifilename);
%     handles = fct_updatedisplay(handles);;
%     %Update handles structure
%     guidata(hObject, handles);
%     
% end
% --------------------------------------------------------------------
function IMPORT_TIF_Callback(hObject, eventdata, handles)
% hObject    handle to IMPORT_TIF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ifilename,ipathname]=uigetfile({'*.tif'},'Film to measure');

if ~strcmp(class(ifilename),'double')
    
    ACTUALDIR=cd;

    channel = fct_colortochannel(handles.defaultcolor);
    button = handles.defaultcolor;
    [z,DPIX,BITS,channel] = fct_read_tif_image(fct_makecleanfilename(ipathname,ifilename),button);

    %     z = fliplr(z);
    if(channel~=-1)
        
        %Initial main variables
        handles.SATMINini = 0;
        handles.SATMAXini = 0;
        handles.BCKGRNDini = [];
        handles.NpixBCKGRNDini = 0;
        handles.DOSEREFini = [];
        handles.NpixDOSEREFini = [];
        handles.DELTAini = DPIX;
        handles.channelini = channel;
        handles.defaultcolor = fct_channeltocolor(channel);
        handles.filterini = 0;
        handles.ORIGINini = [0 0 ]';
        handles.Z = -log10(max(double(z),1)./(2^BITS-1));
%         if channel ==4
%             handles.dataini = 5;
%         else
%             handles.dataini = 2;
%         end
        handles.dataini = 2
        handles.cfilenameini = '';
        handles.BITSini = BITS;
        %Actual main variables
        handles.data = handles.dataini;
        handles.SATMIN = handles.SATMINini;
        handles.SATMAX = handles.SATMAXini;
        %handles.BCKGRND = handles.BCKGRNDini;
        %handles.NpixBCKGRND = handles.NpixBCKGRNDini;
        %handles.DOSEREF = handles.DOSEREFini;
        %handles.NpixDOSEREF = handles.NpixDOSEREFini;
        handles.DELTA = handles.DELTAini;
        handles.channel = handles.channelini;
        handles.filter = handles.filterini;
        handles.ORIGIN = handles.ORIGINini;
        handles.z = handles.Z;
        handles.cfilename = handles.cfilenameini;
        %Other variables
        if handles.channel==4
            handles.color = 'gray';
        else
            handles.color ='jet';
        end
        handles.gridflag = 0;
        handles.gridres = fct_autoruler(handles.z,handles.DELTA);
        handles.ifilename = fct_makecleanfilename(ipathname,ifilename);
    
        handles = fct_updatedisplay(handles);;
        %Update handles structure
        guidata(hObject, handles);
        
    else
        error('Wrong image type.');
    end
    
end

% --------------------------------------------------------------------
function IMPORT_MATLAB_Callback(hObject, eventdata, handles)
% hObject    handle to IMPORT_MATLAB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ifilename,ipathname]=uigetfile({'*.mat'},'Work to import');
if ~strcmp(class(ifilename),'double')
    
    ACTUALDIR = cd;
    filename = fct_makecleanfilename(ipathname,ifilename);
    load (filename);
    if exist('handles_save')
        %I added the following line on 9 March 2015 to fix a bug. The 2014b version of
        %matlab seems to have change the handles on figures
        handles = handles_save     
        close;
        %Initial main variables to former actual variables
        handles.dataini = handles.data;
        handles.SATMINini = handles.SATMIN;
        handles.SATMAXini = handles.SATMAX;
        handles.BCKGRNDini = handles.BCKGRND;
        handles.NpixBCKGRNDini = handles.NpixBCKGRND;
        handles.DOSEREFini = handles.DOSEREF;
        handles.NpixDOSEREFini = handles.NpixDOSEREF;
        handles.DELTAini = handles.DELTA;
        handles.channelini = handles.channel;
        handles.filterini = handles.filter;
        handles.ORIGINini = handles.ORIGIN;
        handles.Z = handles.z;
        handles.cfilenameini = handles.cfilename;
        %Other variables
        if handles.channel==4
            handles.color = 'gray';
        else
            handles.color ='jet';
        end    
        handles.gridflag = 0;
        handles.gridres = fct_autoruler(handles.z,handles.DELTA);
        handles.ifilename = fct_makecleanfilename(ipathname,ifilename);
        handles.Z = handles.z;
        handles.dataini = handles.data;

        H = gcf;
        handles.menus = guihandles(H);
        handles = fct_updatedisplay(handles);
        %Update handles structure
        guidata(hObject, handles);
    end
end

% --------------------------------------------------------------------
function FILTER_RECT_Callback(hObject, eventdata, handles)
% hObject    handle to FILTER_RECT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (fct_isthereanimage(handles)&&(handles.filter==0))
    
    %CHOOSE SIZE OF CONVOLUTION
    width_rect  = inputdlg({'Region width in mm:'},'Rectangular region',1,{'1'});
    width_rect=str2double(width_rect)/10;
    %OPERATIONS
    nx=min(max(ceil(width_rect/handles.DELTA),1),size(handles.z,2));
    ny=min(max(ceil(width_rect/handles.DELTA),1),size(handles.z,1));
    n=min(nx,ny);
    if (nx==size(handles.z,2)||ny==size(handles.z,1))
        error('Filter size greater or equal to image.');
    end
    F=ones(n,n);
    F=F/sum(sum(F));
    %Correction
    z=conv2(handles.z,F,'valid');
    %CREATION OF NEW VARIABLES
    handles.z=z;
    handles.filter=handles.filter+width_rect;
    %UPDATE DISPLAY
    handles = fct_updatedisplay(handles);;
    % Update handles structure
    guidata(hObject, handles);
    
elseif (handles.filter~=0)
    msgbox('You''ve already filtered the image','Warning');
    
end

% --------------------------------------------------------------------
function FILTER_GAUSS_Callback(hObject, eventdata, handles)
% hObject    handle to FILTER_GAUSS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (fct_isthereanimage(handles)&&(handles.filter==0))
    
    %CHOOSE SIZE OF CONVOLUTION
    width_rect  = inputdlg({'Region FWHM (positive) or sigma (negative) in mm:'},'Gaussian',1,{'1'});
    width_rect = str2double(width_rect)/10;
    if width_rect<0
        sigma = -width_rect;
        width_rect = sigma *sqrt(8*log(2));
    else
        sigma = width_rect/sqrt(8*log(2));
    end
    %OPERATIONS
    nx=min(max(ceil(3*sigma/handles.DELTA),1),size(handles.z,2));
    ny=min(max(ceil(3*sigma/handles.DELTA),1),size(handles.z,1));
    n=min(nx,ny);
    if (nx==size(handles.z,2)||ny==size(handles.z,1))
        error('Filter size greater or equal to image.');
    end
    F=ones(n,n);
    for i=1:n
        for j=1:n
            F(i,j)=exp( -(  sqrt( (i-(n-1)/2-1)^2+(j-(n-1)/2-1)^2 ) .* handles.DELTA / sigma  )^2 );
        end
    end
    F=F/sum(sum(F));
    z=conv2(handles.z,F,'valid');
    
    %CREATION OF NEW VARIABLES
    handles.z=z;
    handles.filter = handles.filter + width_rect;
    %UPDATE DISPLAY
    handles = fct_updatedisplay(handles);;
    
    % Update handles structure
    guidata(hObject, handles);
    
elseif (handles.filter~=0)
    msgbox('You''ve already filtered the image','Warning');
end

% --------------------------------------------------------------------
function SHOWSIGNAL_Callback(hObject, eventdata, handles)
% hObject    handle to SHOWSIGNAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS
    z = handles.z;
    h = handles.H;
    
%     if (handles.channel == 4)
%          msgbox('Impossible to obtain signal from multichannel image.','Error');
%          h = gcf;
%     elseif (handles.data == 1)
    if (handles.data == 1)
    elseif (handles.data == 2)
        z = 10.^(-z).*(2^(handles.BITSini)-1);
        handles.data = 1;
    elseif (handles.data == 3)
        msgbox('Impossible to obtain raw signal from net optical density.','Error');
        h = gcf;
    elseif (handles.data == 4)
        msgbox('Impossible to obtain raw signal from dose.','Error');
        h = gcf;
    else
        error('Unsolved error');
    end
    
    %CREATION OF NEW VARIABLES
    handles.z = z;
    
    %UPDATE DISPLAY
    figure(handles.H);
    handles = fct_updatedisplay(handles);;
    figure(h);
    
    % Update handles structure
    guidata(hObject, handles);
    
end

% --------------------------------------------------------------------
function SHOWRAWOD_Callback(hObject, eventdata, handles)
% hObject    handle to SHOWRAWOD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS
    z=handles.z;
    h=handles.H;
    
%     if (handles.channel == 4)
%          msgbox('Impossible to obtain optical density from multichannel image.','Error');
%          h = gcf;
%     elseif (handles.data == 1)
    if (handles.data == 1)
        z = -log10 (max(z,1)./(2^(handles.BITSini)-1)) ;
        handles.data = 2;
    elseif (handles.data == 2)
    elseif (handles.data == 3)
        msgbox('Impossible to obtain raw optical density from net optical density.','Error');
        h = gcf;
    elseif (handles.data == 4)
        msgbox('Impossible to obtain raw optical density from dose.','Error');
        h = gcf;
    else
        error('Unsolved error');
    end
    
    %CREATION OF NEW VARIABLES
    handles.z = z;
    
    %UPDATE DISPLAY
    figure(handles.H);
    handles = fct_updatedisplay(handles);;
    figure(h);
    
    % Update handles structure
    guidata(hObject, handles);
    
end
% --------------------------------------------------------------------
function SHOWNETOD_Callback(hObject, eventdata, handles)
% hObject    handle to SHOWNETOD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
%     if (handles.channel == 4)
%          msgbox('Impossible to obtain net optical density from multichannel image.','Error');
%          h = gcf;
%     elseif handles.data==1
    if handles.data==1
        msgbox('Cannot perform operation from signal.','Error');
    elseif handles.data==3
        msgbox('Operation already performed.','Error');
    elseif handles.data==4
        msgbox('Cannot perform operation from dose.','Error');
    else
        if max(size(handles.BCKGRND))==0
            msgbox('Define background first.','Background');
        else
            %Ici, rien ne prot�ge l'utilisateur d'utiliser un mauvais type
            %de film pour soustraire le background. Par contre, puisque la
            %valeur peut �tre entr�e � la main, la s�curit� est lais�e � la
            %discr�tion de l'usager.
            
            %OPERATIONS & CREATION OF NEW VARIABLES
            if max(size(handles.BCKGRND))==1
                nmb = 1;
            else
                answer  = inputdlg({'Enter the index of the background value to be used:'},'Index',1);
                nmb  = str2double(answer);
            end
            if( (nmb >= 1) && ( nmb <= max(size(handles.BCKGRND)) ) )
                handles.z = max(handles.z-handles.BCKGRND(nmb),0);
                dmb=sprintf('Subtracted value: %.3f',handles.BCKGRND(nmb));
                %dmb=sprintf('Valeur soustraite: %.3f +/- %.3f',handles.BCKGRND(nmb),handles.SBCKGRND(nmb));
                msgbox(dmb,'Background');
                h=gcf;
                handles.data = handles.data + 1;
                
                %UPDATE DISPLAY
                figure(handles.H);
                handles = fct_updatedisplay(handles);;
                figure(h);
                
                % Update handles structure
                guidata(hObject, handles);
            else
                msgbox('Index value not allowed','Error');
            end
        end
    end
end

% --------------------------------------------------------------------
function SHOWDOSE_Callback(hObject, eventdata, handles)
% hObject    handle to SHOWDOSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS
    z = handles.z;
    h = handles.H;
    
    if (handles.data == 1)
        msgbox('Convert to net optical density first.','Error');
        h = gcf;
    elseif (handles.data == 2)
        if max(size(handles.BCKGRND))==0
            msgbox('Define background first.','Bruit de fond');
        else
            %Herem nothing protects the user from using wrong time of film
            %to subtract background. However, since value can be manually
            %entered, security is to the user discretion.
            
            %OPERATIONS & CREATION OF NEW VARIABLES
            if max(size(handles.BCKGRND))==1
                nmb = 1;
            else
                answer  = inputdlg({'Enter the index of the background value to be used:'},'Index',1);
                nmb  = str2double(answer);
            end
            if( (nmb >= 1) && ( nmb <= max(size(handles.BCKGRND)) ) )
                handles.z = max(handles.z-handles.BCKGRND(nmb),0);
                dmb=sprintf('Subtracted value: %.3f',handles.BCKGRND(nmb));
                handles.data = handles.data + 1;
                
                %UPDATE DISPLAY
                figure(handles.H);
                handles = fct_updatedisplay(handles);;
                %figure(h);
                
                % Update handles structure
                guidata(hObject, handles);
            else
                msgbox('Index value not allowed','Error');
            end
        end
        z = handles.z;
        h = handles.H;
        
    elseif (handles.data == 3)
        
        %CREATION OF NEW VARIABLES
        if max(size(handles.cfilename))==0
            
            [ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');
            
            if ~strcmp(class(ifilename),'double')
                file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');
                handles.cfilename = fct_makecleanfilename(ipathname,ifilename);
            end
        else
            file = fopen(handles.cfilename,'r');
        end

        %HB: 30 March 2020: there is a bug if I press cancel to choosing a
        %calibration file
        if file~=-1
            [DOSE,OD,M,opt,sigparam,Npix] = fct_readcalfile(file);
            fclose(file);
            
            [odi,dosei] = fct_getcalcurvepoints(DOSE,OD,M,opt);
            if max(max(z))>max(OD)||min(min(z))<min(OD)
                h_msg = msgbox('Some values of OD are out of range, Values will be saturated','Warning');
                uiwait(h_msg);
                figure(h);
                %z = min(z,max(OD));
                %z = max(z,min(OD));
                z = min(z,max(odi));
                z = max(z,min(odi));
                if handles.SATMAX~=0
                    handles.SATMAX = max(max(z));
                    handles.SATMIN = min(min(z));
                end
            end
            handles.z = interp1(odi,dosei,z);
            if handles.SATMAX~=0
                handles.SATMAX = interp1(odi,dosei,handles.SATMAX);
                handles.SATMIN = interp1(odi,dosei,handles.SATMIN);
            end
            handles.data = 4;
            
            %UPDATE DISPLAY
            figure(handles.H);
            handles = fct_updatedisplay(handles);;
            % Update handles structure
            guidata(hObject, handles);
                       
        end
    elseif (handles.data == 5)
        %CREATION OF NEW VARIABLES
        if max(size(handles.cfilename))==0
            
            [ifilename,ipathname] = uigetfile({'*.mlt'},'Choose calibration curve');
            
            if ~strcmp(class(ifilename),'double')
                file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');
                handles.cfilename = fct_makecleanfilename(ipathname,ifilename);
            else 
                file = -1;
            end
        else
            file = fopen(handles.cfilename,'r');
        end
        if file~=-1
            [DOSE,XI,sXI,c,V,rescal,sximesh] = fct_ReadCalFileMulti(file);
            fclose(file);
            
            handles.z = (z-c(1))/c(2);
            handles.data = 4;
            
            %UPDATE DISPLAY
            figure(handles.H);
            handles = fct_updatedisplay(handles);;
            % Update handles structure
            guidata(hObject, handles);
                       
        end        
    elseif handles.data==4
    else
        error('Unsolved error');
    end    
end


% --------------------------------------------------------------------
function COLORMAP_Callback(hObject, eventdata, handles)
% hObject    handle to COLORMAP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    
    str{1} = 'winter';
    str{2}='gray';
    str{3}='autumn';
    str{4}='summer';
    str{5}='spring';
    str{6}='bone';
    str{7}='colorcube';
    str{8}='cool';
    str{9}='copper';
    str{10}='flag';
    str{11}='hot';
    str{12}='hsv';
    str{13}='lines';
    str{14}='pink';
    str{15}='prism';
    str{16}='jet';
    
    [clmp,ok] = listdlg('Name','Colormap','ListString',str);
    
    if (ok==1)&&(length(clmp)==1)
        
        clmp = str{clmp};
        %CREATION OF NEW VARIABLES
        handles.color = clmp;
        %UPDATE DISPLAY
        handles = fct_updatedisplay(handles);;
        % Update handles structure
        guidata(hObject, handles);
        
    end
    
end

% --------------------------------------------------------------------
function GRIDSIZE_Callback(hObject, eventdata, handles)
% hObject    handle to gridSIZE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  
if fct_isthereanimage(handles)
    
    v = [1 2 5 10 20 50];
    
    for i=1:length(v)
        str{i} = sprintf('%d',v(i));
    end
    [res,ok] = listdlg('Name','Grid size (mm)','ListString',str);
    
    if (ok==1)&&(length(res)==1)
        
        res = str2num(str{res})/10;
        handles.gridres = res;
        %UPDATE DISPLAY
        if fct_isthereanimage(handles)
            handles = fct_updatedisplay(handles);;
        end
        % Update handles structure
        guidata(hObject, handles);
    end
end


% --------------------------------------------------------------------
function VERSION_Callback(hObject, eventdata, handles)
% hObject    handle to VERSION (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_version()==1
    close all;
end

% --------------------------------------------------------------------
function CALCURVE_CREATE_Callback(hObject, eventdata, handles)
% hObject    handle to CALFILM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% button = questdlg('What type of calibration do you wish?','Channel','Single','Multi','Multi');
button = 'Single';
if strcmp(button,'Multi')
    err = fct_CreateCalCurveMulti(handles);
else
%     err = fct_CreateCalCurveSingle(handles);
    err = fct_CreateCalCurveSingleNewUncertROI(handles);
end

% --------------------------------------------------------------------
function BCKGRND_DEFINE_Callback(hObject, eventdata, handles)
% hObject    handle to BCKGRND_DEFINE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)

%     if (handles.channel == 4)
%          msgbox('Cannot perform operation from multichannel image.','Error');
%          h = gcf;
%     elseif handles.data==1
    if handles.data==1
        msgbox('Cannot perform operation from signal.','Signal');
    elseif handles.data==3
        msgbox('Cannot perform operation from net optical density.','OD nette');
    elseif handles.data==4
        msgbox('Cannot perform operation from dose.','Dose');
    else
        %OPERATIONS & CREATION OF NEW VARIABLES
        answer  = inputdlg({'Enter number of background values to be defined:'},'Number of values',1);
        nmb  = str2double(answer);
        if max(size(answer))==0
        else
            nmb  = str2double(answer);
            handles.BCKGRND = zeros(nmb,1);
            button = questdlg('How do you wish to define background values?','Background definition','Region','Manually','Region');
            if strcmp(button,'Manually')
                for i=1:nmb
                    prompt{1} = 'Enter background value:';
                    %prompt{2} = 'Enter number of pixels';
                    answer  = inputdlg(prompt,'Background',1);
                    answer = str2double(answer);
                    handles.BCKGRND(i)  = answer(1);
                    handles.NpixBCKGRND = answer(1)*0 + 1000000000;
                end
            else
                if 1
                    button = questdlg('What type of calibration do you wish to perform','Type','Fast','Precise','Precise') ;
                else
                    button = 'Precise';
                end
                if strcmp(button,'Precise')
                    short = 0;
                else
                    short = 1;
                end
                width_rect  = inputdlg({'Region width in cm:','Region height in cm:'},'Rectangular region',1,{'1','1'})
                deltax = str2double(width_rect(1));
                deltay = str2double(width_rect(2));
                handles.NpixBCKGRND = deltax*deltay/handles.CCDres/handles.CCDres;
                if short
                    handles.BCKGRND = zeros(nmb,1);
                    dumb = sprintf('Click at the center of each of the %d region of interest in order and press enter after each point.',nmb);
                    msg = questdlg(dumb,'Unirradiated films','OK','OK');
                    figure(handles.H);
                    for i=1:nmb
                        [newz,center,owidth] = fct_getroi(handles.z,handles.DELTA,'Point',[deltax deltay]);
                        [mu, sigmamu, s] = fct_analyze_region(newz);
                        handles.BCKGRND(i) = mu
                    end
                else
                    figure(handles.H);
                    for i=1:nmb
                        [newz,center,owidth] = fct_getroi(handles.z,handles.DELTA,'Zoom-fixed',[deltax deltay]);
                        [mu, sigmamu, s] = fct_analyze_region(newz);
                        handles.BCKGRND(i) = mu
                    end
                end
                if nmb>1
                    button = questdlg('What kind of background analysis do you want to do?','Background analysis','Average','One-by-one','One-by-one') ;
                    if strcmp(button,'Average')
                        %ici ajouter un while au cas ou cela ne marche pas
                        handles.BCKGRND = mean(handles.BCKGRND);
                    end
                end
            end
        end
        h=gcf;
        %UPDATE DISPLAY
        figure(handles.H);
        handles = fct_updatedisplay(handles);;
        figure(h);
        
        %HB 15 April 2015: this is the new way of dealing with background
        %we take the average during calibration and we assume the ROI will be
        %large enough during use
        handles.NpixBCKGRND = 1e10;
        
        % Update handles structure
        guidata(hObject, handles);
    end
end

% --------------------------------------------------------------------
function DEFINE_REF_DOSE_Callback(hObject, eventdata, handles)
% hObject    handle to DEFINE_REF_DOSE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    if handles.data==1
        msgbox('Cannot perform operation from signal.','Signal');
    elseif handles.data==2
        msgbox('Cannot perform operation from raw optical density.','Raw OD');
    elseif handles.data==3
        msgbox('Cannot perform operation from net optical density.','Net OD');
    elseif handles.data==5
        msgbox('Cannot perform operation from xi.','Xi');
    else
        button = questdlg('How do you wish to define the value?','Reference dose','Region','Manually','Region');
        if strcmp(button,'Region')
            [A,center,owidth] = fct_getroi(handles.z,handles.DELTA,'Rectangular',0);
            [mu, sigma, Npix] = fct_analyze_region(A);
            Npix = Npix*(handles.DELTA/handles.CCDres)^2;
            if handles.channel==4
                file = fopen(handles.cfilename,'r');
                [DOSE,XI,sXI,c,V,rescal,sximesh] = fct_ReadCalFileMulti(file);
                fclose(file);
                [sd,sxi] = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,Npix,mu);
                dmb = sprintf('Reference DI = %.1f +/- %.1f (type A, k=1)\nROI has %.0f elem. pixels',mu,sd,Npix);
                hinfo = helpdlg(dmb,'Reference');
                figure(hinfo);
                uiwait(hinfo);
                handles.DOSEREF  = mu;
                handles.NpixDOSEREF = Npix;
            else
                V = fct_getcovarmatrix(mu,Npix,handles.NpixBCKGRND,handles.cfilename);
                smu = sqrt(V);
                dmb = sprintf('Reference value is %.1f +/- %.1f (%.1f%%)\nSelected region has %.0f pixels',mu,smu,smu/mu*100,Npix);
                hinfo = helpdlg(dmb,'Reference');
                figure(hinfo);
                uiwait(hinfo);
                handles.DOSEREF  = mu;
                handles.NpixDOSEREF = Npix;
            end
        else
            prompt{1} = 'Enter DI value:';
            prompt{2} = 'Enter nb. of elem. pixels:';
            answer  = inputdlg(prompt,'Reference dose index',1,{'200','1e09'});
            answer
            if numel(answer)~=2
            else
                answer = str2double(answer);
                if (answer(1)<0)||(answer(2)<0)||(sum(isnan(answer))~=0) 
                    errordlg('Invalid entry. Values not defined.');
                else
                    handles.DOSEREF  = answer(1);
                    handles.NpixDOSEREF  = floor(answer(2));
                end
            end
        end
        h=gcf;
        %UPDATE DISPLAY
        figure(handles.H);
        handles = fct_updatedisplay(handles);;
        figure(h);
        
        % Update handles structure
        guidata(hObject, handles);
        
    end
end


% --------------------------------------------------------------------
function SATURATE_MAX_Callback(hObject, eventdata, handles)
% hObject    handle to SATURATE_MAX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    if handles.data==1
        msgbox('Cannot perform operation from signal.','Signal');
    else
        %OPERATIONS & CREATION OF NEW VARIABLES
        button = questdlg('How do you wish to define the saturation?','Definition of saturation','Region','Manually','Region');
        if strcmp(button,'Region')
            figure(handles.H);
            %A=imcrop;
            A = fct_getroi(handles.z,handles.DELTA,'Free',0);
            handles.SATMAX = max(max(A));
        else
            answer  = inputdlg({'Enter saturation value:'},'Saturation',1);
            handles.SATMAX = str2double(answer);
        end
        
        h = gcf;
        handles.z = min(handles.z,handles.SATMAX);
        %guidata(hObject, handles);
        
        %UPDATE DISPLAY
        figure(handles.H);
        handles = fct_updatedisplay(handles);;
        figure(h);
        
        % Update handles structure
        guidata(hObject, handles);
    end
end

% --------------------------------------------------------------------
function SATURATE_MIN_Callback(hObject, eventdata, handles)
% hObject    handle to SATURATE_MIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    if handles.data==1
        msgbox('Cannot perform operation from signal.','Signal');
    else
        %OPERATIONS & CREATION OF NEW VARIABLES
        button = questdlg('How do you wish to define saturation?','Definition of saturation','Region','Manually','Region');
        if strcmp(button,'Region')
            figure(handles.H);
            %A=imcrop;
            A = fct_getroi(handles.z,handles.DELTA,'Free',0);
            handles.SATMIN = min(min(A));
        else
            answer  = inputdlg({'Enter saturation value:'},'Saturation',1);
            handles.SATMIN = max(0,str2double(answer));
        end
        if handles.SATMAX==0
            handles.SATMAX = max(max(handles.z));
        end
        h=gcf;
        handles.z = max(handles.z,handles.SATMIN);
        
        %UPDATE DISPLAY
        figure(handles.H);
        handles = fct_updatedisplay(handles);;
        figure(h);
        
        % Update handles structure
        guidata(hObject, handles);
    end
end

% --------------------------------------------------------------------
function MEANDOSE_ORIGIN_Callback(hObject, eventdata, handles)
% hObject    handle to MEANDOSE_ORIGIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    answer  = inputdlg({'Region width in mm:','Region height in mm:', ...
        'Maximum error between measured origin and real position in mm'},'Rectangular region',1,{'1','1','1'});
    if numel(answer)==0
    else
        answer = str2double(answer)/10;
        width_rect = [answer(1) answer(2)]';
        maxshift = answer(3);

        nlines = size(handles.z,1);
        ncols = size(handles.z,2);
        [x,y] = fct_gridindextopos(nlines,ncols,handles.DELTA);
        x = x - handles.ORIGIN(1);
        y = y - handles.ORIGIN(2);

        %get value at believed origin
        k1 = intersect(find(x > -width_rect(1)/2),find(x < width_rect(1)/2) );
        k2 = intersect(find(y > -width_rect(2)/2),find(y < width_rect(2)/2) );
        A = handles.z(k2,k1);
        [mu, sigma, Npix] = fct_analyze_region(A);
        Npix = Npix*(handles.DELTA/handles.CCDres)^2;
        [sizeROI1,sizeROI2] = size(A);
        sizeROI1 = sizeROI1*handles.DELTA;
        sizeROI2 = sizeROI2*handles.DELTA;
        exp = ['^2'];
        dumb1=sprintf('ROI = %.2f x %.2f mm%s',sizeROI2*10,sizeROI1*10,exp);
        %estimate variance shifting origin with r normally distributed
        %and phi equally distributed
        if maxshift==0
            sigma(1:3) = 0;
        else
            sigma(1:3) = 0;
            K = 500;
            s(1:3) = 0;
            s2(1:3) = 0;
            for i=1:K
                r = 3/2*maxshift;
                while r>maxshift
                    r = fct_randgauss(0,maxshift/2,1);
                end
                phi = 2*pi*rand(1,1);
                x0 = r*cos(phi);
                y0 = r*sin(phi);
                k1 = intersect(find(x > x0-width_rect(1)/2),find (x < x0+width_rect(1)/2) );
                k2 = intersect(find(y > y0-width_rect(2)/2),find (y < y0+width_rect(2)/2) );
                A = handles.z(k2,k1);
                [mui, sigmai, Ni] = fct_analyze_region(A);
                s(1)  = s(1)  + mui;
                s2(1) = s2(1) + mui*mui;

                r = fct_randrect(0,maxshift,1);
                phi = 2*pi*rand(1,1);
                x0 = r*cos(phi);
                y0 = r*sin(phi);
                k1 = intersect(find(x > x0-width_rect(1)/2),find (x < x0+width_rect(1)/2) );
                k2 = intersect(find(y > y0-width_rect(2)/2),find (y < y0+width_rect(2)/2) );
                A = handles.z(k2,k1);
                [mui, sigmai, Ni] = fct_analyze_region(A);
                s(2)  = s(2)  + mui;
                s2(2) = s2(2) + mui*mui;

                r = fct_randtriang(0,maxshift,1);
                phi = 2*pi*rand(1,1);
                x0 = r*cos(phi);
                y0 = r*sin(phi);
                k1 = intersect(find(x > x0-width_rect(1)/2),find (x < x0+width_rect(1)/2) );
                k2 = intersect(find(y > y0-width_rect(2)/2),find (y < y0+width_rect(2)/2) );
                A = handles.z(k2,k1);
                [mui, sigmai, Ni] = fct_analyze_region(A);
                s(3)  = s(3)  + mui;
                s2(3) = s2(3) + mui*mui;
            end
            sigma = sqrt(s2/(K-1)-s.*s/K/(K-1));
        end
        
        if (handles.data==1)
            dumb = sprintf('%s\n\nValue = %.0f\n\nUncertainty from positionning:\nGaussian = +/- %.0f\nCircular = +/- %.0f\nTriangular = +/- %.0f',dumb1,mu,sigma(1),sigma(2),sigma(3));
        elseif (handles.data==2)||(handles.data==3)||(handles.data==5)
            dumb = sprintf('%s\n\nValue = %.3f\n\nUncertainty from positionning:\nGaussian= +/- %.3f\nCircular = +/- %.3f\nTriangular = +/- %.3f',dumb1,mu,sigma(1),sigma(2),sigma(3));
        elseif (handles.data==4)
            if (handles.channel==4)
                file = fopen(handles.cfilename,'r');
                [DOSE,XI,sXI,c,V,rescal,sximesh] = fct_ReadCalFileMulti(file);
                fclose(file);
                if  max(size(handles.DOSEREF))==0;
                    [sd,sxi] = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,Npix,mu);
                    dumb=sprintf('%s\n\nDI = %.1f +/- %.1f (type A, k=1)',dumb1,mu,sd);
                else
                    [sd,sxi] = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,Npix,mu);
                    dumb = sprintf('%s\n\nDI = %.1f +/- %.1f (type A, k=1)',dumb1,mu,sd);
                    r = mu/handles.DOSEREF;
                    sr = fct_UncertaintyMultichannelSingleMeasRel(c,V,sximesh,handles.NpixDOSEREF,handles.DOSEREF,Npix,mu);
                    dumb = sprintf('%s\nDI ratio = %.3f +/- %.3f (type A, k=1)',dumb,r,sr);
                end            
                dumb = sprintf('%s\n\nUncertainty from positionning:\nGaussian = +/- %.1f (%.1f %%)\nCircular = +/- %.0f (%.1f %%)\nTriangular = +/- %.0f (%.1f %%)',...
                    dumb,sigma(1),sigma(1)/mu*100,sigma(2),sigma(2)/mu*100,sigma(3),sigma(3)/mu*100);
            else
                if max(size(handles.DOSEREF))==0;
                    DOSE = mu;
                else
                    DOSE =[mu handles.DOSEREF];
                    Npix = [Npix handles.NpixDOSEREF];
                end
                V = fct_getcovarmatrix(DOSE,Npix,handles.NpixBCKGRND,handles.cfilename);
                dumb = fct_meanval_output(DOSE,V);
                dumb = sprintf('%s\n\nUncertainty from positionning:\nGaussian = +/- %.1f (%.1f %%)\nCircular = +/- %.0f (%.1f %%)\nTriangular = +/- %.0f (%.1f %%)',...
                    dumb,sigma(1),sigma(1)/mu*100,sigma(2),sigma(2)/mu*100,sigma(3),sigma(3)/mu*100);
            end
        end
        msgbox(dumb,'Mean');
        h=gcf;
        figure(handles.H);
        handles = fct_updatedisplay(handles);;
        figure(h);
    end
end

% --------------------------------------------------------------------
function MEANDOSE_SELECT_Callback(hObject, eventdata, handles)
% hObject    handle to MENU_MEANDOSE_SELECT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)    
    
    str{1} = 'Free';
    str{2} = 'Rectangular';
    str{3} = 'Zoom';
    str{4} = 'Zoom-Fixed';
    str{5} = 'Fixed';
    str{6} = 'Point';
    [type,ok] = listdlg('Name','ROI','ListString',str);

    if (ok==1)&&(length(type)==1)
        type = str{type};
        figure(handles.H);
        [A,center,owidth] = fct_getroi(handles.z,handles.DELTA,type,0);
        %OPERATIONS
        [mu, sigma, Npix] = fct_analyze_region(A);
        if ~isnan(mu)            
            Npix = Npix*(handles.DELTA/handles.CCDres)^2;
            [sizeROI1,sizeROI2] = size(A);
            sizeROI1 = sizeROI1*handles.DELTA;
            sizeROI2 = sizeROI2*handles.DELTA;
            exp = ['^2'];
            dumb1=sprintf('ROI = %.2f x %.2f mm%s',sizeROI2*10,sizeROI1*10,exp);
            if (handles.data==1)
                dumb=sprintf('%s\n\nMean = %.0f\nStandard deviation = %.0f\nNb. elem. pixels = %d',dumb1,mu,sigma,Npix);
            elseif (handles.data==2)||(handles.data==3)||(handles.data==5)
                dumb=sprintf('%s\n\nMean = %.4f\nStandard deviation = %.4f\nNb. elem. pixels = %d',dumb1,mu,sigma,Npix);
            elseif (handles.data==4)
                if (handles.channel==4)
                    file = fopen(handles.cfilename,'r');
                    [DOSE,XI,sXI,c,V,rescal,sximesh] = fct_ReadCalFileMulti(file);
                    fclose(file);
                    if  max(size(handles.DOSEREF))==0;
                        sd = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,Npix,mu);
                        dumb=sprintf('%s\n\nDI = %.1f +/- %.1f (type A, k=1)',dumb1,mu,sd);
                    else
                        sd = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,Npix,mu);
                        dumb = sprintf('%s\n\nDI = %.1f +/- %.1f (type A, k=1)',dumb1,mu,sd);
                        r = mu/handles.DOSEREF;
                        sr = fct_UncertaintyMultichannelSingleMeasRel(c,V,sximesh,handles.NpixDOSEREF,handles.DOSEREF,Npix,mu);
                        dumb = sprintf('%s\nDI ratio = %.3f +/- %.3f (type A, k=1)',dumb,mu/handles.DOSEREF,sr);
                    end
                else
                    if max(size(handles.DOSEREF))==0;
                        DOSE = mu;
                    else
                        DOSE =[mu handles.DOSEREF];
                        Npix = [Npix handles.NpixDOSEREF];
                    end
                    V = fct_getcovarmatrix(DOSE,Npix,handles.NpixBCKGRND,handles.cfilename);
                    dumb = fct_meanval_output(DOSE,V);                
                end
            else
                errordlg('Unsolved error.');
            end
            msgbox(dumb,'Mean');
            h=gcf;

            %UPDATE DISPLAY
            figure(handles.H);
            handles = fct_updatedisplay(handles);;
            figure(h);

            % Update handles structure
            %guidata(hObject, handles);
        end
    end
end

% --------------------------------------------------------------------
function MEANDOSE_MULTIPLE_Callback(hObject, eventdata, handles)
% hObject    handle to MEANDOSE_POINT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if fct_isthereanimage(handles)
    
    h = figure(handles.H);
    %fig_h = figure;
    %fct_updatedisplay(handles);
    
    str{1}='Rectangular';
    str{2}='Zoom';
    str{3}='Zoom-Fixed';
    str{4}='Fixed';
    str{5}='Point';
    
    [type,ok] = listdlg('Name','ROI','ListString',str);
    
    if (ok==1)
        type = str{type};
        
        nbfilms  = inputdlg({'Number films'},'Number of films',1);
        if numel(nbfilms)==0
        else
            nbfilms  = str2double(nbfilms);
            width_rect  = inputdlg({'Region width in cm:','Region height in cm:'},'Rectangular region',1,{'1','1'});
            if numel(width_rect)==0
            else
                iwidth = [str2double(width_rect(1)) str2double(width_rect(2))];
                val(1:nbfilms) = 0;

                flag = 0;
                i = 0;
                while i <nbfilms
                    i = i+1;
                    [A,center,owidth] = fct_getroi(handles.z,handles.DELTA,type,iwidth);
                    [m1,m2] = size(A);
                    [sizeROI1,sizeROI2] = size(A);
                    sizeROI1 = sizeROI1*handles.DELTA;
                    sizeROI2 = sizeROI2*handles.DELTA;
                    exp = ['^2'];
                    dumb1=sprintf('ROI = %.2f x %.2f mm%s',sizeROI2*10,sizeROI1*10,exp);
                    if m1*m2==0
                        flag = 1
                        i = nbfilms;
                    else
                        [mu, sigma, Npix] = fct_analyze_region(A);
                        val(i) = mu;
                        Npix = Npix*(handles.DELTA/handles.CCDres)^2;
                        if (handles.data==1)
                            dumb{i} = sprintf('Value #%d\n%s\nMean = %.0f\nStandard deviation = %.0f\nNb. elem. pixels = %d',i, dumb1, mu,sigma,Npix);
                        elseif (handles.data==2)||(handles.data==3)||(handles.data==5)
                            dumb{i} = sprintf('Value #%d\n%s\nMean = %.4f\nStandard deviation = %.4f\nNb. elem. pixels = %d',i, dumb1, mu,sigma,Npix);
                        elseif (handles.data==4)
                            if (handles.channel==4)
                                file = fopen(handles.cfilename,'r');
                                [DOSE,XI,sXI,c,V,rescal,sximesh] = fct_ReadCalFileMulti(file);
                                fclose(file);
                           
                                if  max(size(handles.DOSEREF))==0;
                                    sd = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,Npix,mu); 
                                    dumb{i} = sprintf('Value #%d\n%s\nDI = %.1f +/- %.1f (type A, k=1)',i,dumb1,mu,sd);
                                else
                                    sd = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,Npix,mu); 
                                    dumb{i} = sprintf('Value #%d\n%s\nDI = %.1f +/- %.1f (type A, k=1)',i,dumb1,mu,sd);
                                    r = mu/handles.DOSEREF;
                                    sr = fct_UncertaintyMultichannelSingleMeasRel(c,V,sximesh,handles.NpixDOSEREF,handles.DOSEREF,Npix,mu);
                                    str = '_{ref}';
                                    dumb{i} = sprintf('%s\nDI ratio = %.3f +/- %.3f (type A, k=1)',dumb{i},mu/handles.DOSEREF,sr);                                
                                end
                            else        
                                if max(size(handles.DOSEREF))==0;
                                    DOSE = mu;
                                else
                                    DOSE =[mu handles.DOSEREF];
                                    Npix = [Npix handles.NpixDOSEREF];
                                end
                                V = fct_getcovarmatrix(DOSE,Npix,handles.NpixBCKGRND,handles.cfilename);
                                dumb{i} = fct_meanval_output(DOSE,V);
                            end
                        end
                    end
                end
                if flag == 0
                    dmb = sprintf('%s',dumb{1});
                    for i = 2:nbfilms
                        dmb = sprintf('%s\n\n%s',dmb,dumb{i});
                    end
                    val(:);
                    msgbox(dmb,'Mean');
                    h = gcf;
                    %UPDATE DISPLAY
                    figure(handles.H);
                    handles = fct_updatedisplay(handles);;
                    figure(h);
                    uiwait(h);
                    % Update handles structure
                    %guidata(hObject, handles);
                else
                    errordlg('Invalid entry. Values not defined.');
                end
            end
        end
    end
end


% --------------------------------------------------------------------
function SAVE_MATLAB_Callback(hObject, eventdata, handles)
% hObject    handle to SAVE_MATLAB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    tmp = handles.ifilename;
    k = max(find(tmp=='.'));
    if size(k>0)
        tmp = [tmp(1:k-1) '_worked.mat'];
    end
    [ifilename,ipathname]=uiputfile({'*.mat'},'Work to save',tmp);
    ACTUALDIR = cd;
    if ifilename==0
    else
        filename = fct_makecleanfilename(ipathname,ifilename);
        handles_save = handles;
        save(filename,'handles_save');
    end
end

% --------------------------------------------------------------------
function SET_ORIGIN_Callback(hObject, eventdata, handles)
% hObject    handle to SET_ORIGIN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS & CREATION OF NEW VARIABLES
    %button = questdlg('Comment voulez-vous definir l''origine?','Definition de l''origine','Point','Manually','Point');
    button = 'Point';
    if strcmp(button,'Point')
        figure(handles.H);
        [x,y,nx,ny]=fct_getpoints(handles.z,handles.DELTA);
        answer = [x(1) y(1)]';
        handles.ORIGIN = answer +  handles.ORIGIN;
        
    else
        prompt{1} = 'Enter origin value in x:';
        prompt{2} = 'Enter origin value in y:';
        answer  = inputdlg(prompt,'Origin',1);
        handles.ORIGIN  = str2double(answer);
    end
    
    h=gcf;
    %UPDATE DISPLAY
    figure(handles.H);
    handles = fct_updatedisplay(handles);;
    figure(h);
    
    % Update handles structure
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function REGISTRATION_CROSSHAIR_Callback(hObject, eventdata, handles)
% hObject    handle to REGISTRATION_CROSSHAIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    %OPERATIONS & CREATION OF NEW VARIABLES
    figure(handles.H);
    [x,y,nx,ny]= fct_getpoints(handles.z,handles.DELTA);
    if length(x)~=4
        errordlg('Invalid entry. Try 4 points.');
    else
        x = sort(x(1:4));
        y = sort(y(1:4));
        answer = [(x(2) + x(3))/2 (y(2) + y(3))/2]';
        handles.ORIGIN = answer +  handles.ORIGIN;
    %     handles.gridflag = 1;
    %     handles.gridres = 0.1;
        h=gcf;
        %UPDATE DISPLAY
        figure(handles.H);
        handles = fct_updatedisplay(handles);;
        figure(h);

        % Update handles structure
        guidata(hObject, handles);
    end
end

% % --------------------------------------------------------------------
% function REGISTRATION_MANUAL_Callback(hObject, eventdata, handles)
% % hObject    handle to REGISTRATION_MANUAL (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% if fct_isthereanimage(handles)
%     
%     [ifilename,ipathname]=uigetfile({'*.mat'},'Reference fluence map');
%     
%     if (~strcmp(class(ifilename),'double'))
%         
%         ACTUALDIR=cd;
%         filename = fct_makecleanfilename(ipathname,ifilename)
%         load (filename);
%         
%         if exist('X')
%             
%             answer  = inputdlg({'Enter fluence map magnification factor'},'Magnification',1);
%             scale = str2double(answer);
%             X = X*scale;
%             Y = Y*scale;
%             
%             handles.gridflag = 0;
%             handles = fct_updatedisplay(handles);;
%             
%             SATMIN = min(min(handles.z))
%             SATMAX = max(max(handles.z))
%             delta = abs(X(1)-X(2));
%             
%             minZ = min(min(Z))
%             maxZ = max(max(Z))
%             
%             m = (SATMAX-SATMIN)/(maxZ-minZ)
%             b = SATMIN
%             Z = m*(Z-minZ) +b;
%             
%             h = figure;
%             imagesc(X,Y,Z);
%             gridres = 0.1;
%             vect1 = [gridres*ceil(min(X)/gridres)-gridres:gridres:gridres*floor(max(X)/gridres)+gridres];
%             vect2 = [gridres*ceil(min(Y)/gridres)-gridres:gridres:gridres*floor(max(Y)/gridres)+gridres];
%             set(gca,'DataAspectRatio',[1 1 1],'Xtick',vect1,'Ytick',vect2);
%             set(gcf,'Units','pixels');
%             impixelinfo;
%             colormap(handles.color);
%             colorbar('vert');
%             grid off;
%             set(gca,'DataAspectRatio',[1 1 1]);
%             
%             button = questdlg('Do you wish to flip the image horizontally?','Image flip','Yes','No','No');
%             if strcmp(button,'Yes')
%                 Z = fliplr(Z);
%                 imagesc(X,Y,Z);
%                 gridres = 0.1;
%                 vect1 = [gridres*ceil(min(X)/gridres)-gridres:gridres:gridres*floor(max(X)/gridres)+gridres];
%                 vect2 = [gridres*ceil(min(Y)/gridres)-gridres:gridres:gridres*floor(max(Y)/gridres)+gridres];
%                 set(gca,'DataAspectRatio',[1 1 1],'Xtick',vect1,'Ytick',vect2);
%                 set(gcf,'Units','pixels');
%                 impixelinfo;
%                 colormap(handles.color);
%                 colorbar('vert');
%                 grid off;
%                 set(gca,'DataAspectRatio',[1 1 1]);
%             end
%             
%             button = questdlg('Do you wish to flip the image vertically?','Image flip','Yes','No','No');
%             if strcmp(button,'Yes')
%                 Z = flipud(Z);
%                 imagesc(X,Y,Z);
%                 gridres = 0.1;
%                 vect1 = [gridres*ceil(min(X)/gridres)-gridres:gridres:gridres*floor(max(X)/gridres)+gridres];
%                 vect2 = [gridres*ceil(min(Y)/gridres)-gridres:gridres:gridres*floor(max(Y)/gridres)+gridres];
%                 set(gca,'DataAspectRatio',[1 1 1],'Xtick',vect1,'Ytick',vect2);
%                 set(gcf,'Units','pixels');
%                 impixelinfo;
%                 colormap(handles.color);
%                 colorbar('vert');
%                 grid off;
%                 set(gca,'DataAspectRatio',[1 1 1]);
%             end
%             
%             [x,y,nx,ny]=fct_getpoints(flipud(Z),delta);
%             answer1 = [x(1) y(1)]';
%             close(h);
%             figure(handles.H);
%             [x,y,nx,ny]=fct_getpoints(handles.z,handles.DELTA);
%             answer2 = [x(1) y(1)]';
%             handles.ORIGIN = handles.ORIGIN + answer2 - answer1;
%             handles.gridflag = 1;
%             handles.gridres = 0.1;
%             %UPDATE DISPLAY
%             figure(handles.H);
%             handles = fct_updatedisplay(handles);;
%             % Update handles structure
%             guidata(hObject, handles);
%         else
%             msgbox('Fichier non valide','Error');
%             h=gcf;
%         end
%         
%     end
% end

% --------------------------------------------------------------------
function CALCURVE_VISUALIZE_Callback(hObject, eventdata, handles)
% hObject    handle to CALCURVE_CREATE_VIEW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if strcmp(handles.defaultcolor,'Multi')
%     [ifilename,ipathname] = uigetfile({'*.mlt'},'Choose calibration curve');
% else
%     [ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');
% end
[ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');
if ~strcmp(class(ifilename),'double')
    file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');
%     if strcmp(handles.defaultcolor,'Multi')
%         [DOSE,XI,sXI,c,V,rescal,sximesh] = fct_ReadCalFileMulti(file);
%         h = fct_show_calcurve_multi(DOSE,XI,sXI,c,V,rescal,sximesh);
%     else
        [DOSE,OD,M,opt,sigparam,Npix] = fct_readcalfile(file);
        if 1
            button = questdlg('Do you wish to see all functions?','Vizualisation','Yes','No','No');
        else
            button = 'No';
        end
        if strcmp(button,'No')
            show_forcing = 0;
        else
            show_forcing = 1;
        end
        %show_forcing = 0;
        h = fct_show_calcurve(DOSE,OD,M,opt,show_forcing);
%     end
    title(fct_addbackslash(fct_makecleanfilename(ipathname,ifilename)));
    fclose(file);
end

% --------------------------------------------------------------------
function CALCURVE_SETDEFAULT_Callback(hObject, eventdata, handles)
% hObject    handle to CALCURVE_SETDEFAULT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if strcmp(handles.defaultcolor,'Multi')
    [ifilename,ipathname] = uigetfile({'*.mlt'},'Choose calibration curve');
else
    [ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');
end
if ~strcmp(class(ifilename),'double')
    file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');

    if strcmp(handles.defaultcolor,'Multi')
        [DOSE,XI,sXI,c,V,rescal,sximesh] = fct_ReadCalFileMulti(file);
    else
        [DOSE,OD,M,opt,sigparam,Npix] = fct_readcalfile(file);
    end
    fclose(file);
    handles.cfilename = fct_makecleanfilename(ipathname,ifilename);
    
    if fct_isthereanimage(handles)
        %UPDATE DISPLAY
        figure(handles.H);
        handles = fct_updatedisplay(handles);;
    end
    %h = fct_show_calcurve(DOSE,OD,N,opt,0);
    % Update handles structure
    guidata(hObject, handles);
    
end

% --------------------------------------------------------------------
function MEANTIF16RGB_Callback(hObject, eventdata, handles)
% hObject    handle to MEANTIF16RGB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Multiifilename,ipathname] = uigetfile({'*.tif'},'Choose .tif file to be averaged','MultiSelect','on');
S = size(Multiifilename,2);

if (size(char(Multiifilename(:,1)),2) == 1)
    if Multiifilename==0
    else
        msgbox('Must select at least 2 images');
    end
else
    iscell(Multiifilename(:,1))
    ifilename = char(Multiifilename(:,1));
    tfilename = ifilename;
    ifilename = fct_makecleanfilename(ipathname, ifilename);
    h = waitbar(1/S,'Please wait');
    [dummy,Resolution0,bits,channel] = fct_read_tif_image(ifilename,'All');
    IMAGE = double(dummy);
    N = 1;
    errflag = 0;
    i=1;
    while(i<=S)&&(~errflag)
        ifilename = char(Multiifilename(:,i))
        ifilename = fct_makecleanfilename(ipathname,ifilename);
        waitbar(i/S,h);
        [dummy,Resolution,bits,channel] = fct_read_tif_image(ifilename,'All');
        channel
        Resolution
        if (Resolution==0)||(channel~=5)
            close(h);
            errordlg('Image should be 16 bits RGB only.');
            errflag = 1;
        elseif (Resolution~=0)&&(Resolution == Resolution0)&&(size(IMAGE,1)==size(dummy,1))&&(size(IMAGE,2)==size(dummy,2))
            IMAGE = IMAGE + double(dummy);
            N = N +1;
        else
            close(h);
            errordlg([ifilename ' has a different resolution or different size from previous one(s).'],'Error');
            errflag = 1;
        end
        i = i+1;
    end
    if (~errflag)
        close(h);
        IMAGE = double(IMAGE) / N;
        %rotate or flip image
        flg = 0;%this is to skip the following while loop
        while flg
            figure('NumberTitle','off','Name','Averaged image');
            h = fct_display(uint16(IMAGE),round(2.54/Resolution));
            ans = questdlg('Do you want to rotate the image?','Image','90 left','90 right','No','No') ;
            if strcmp(ans,'90 right')
                tIMAGE(:,:,1) = fliplr(IMAGE(:,:,1)');
                tIMAGE(:,:,2) = fliplr(IMAGE(:,:,2)');
                tIMAGE(:,:,3) = fliplr(IMAGE(:,:,3)');
                IMAGE = tIMAGE; clear tIMAGE;
                close(h);figure('NumberTitle','off','Name','Averaged image');
                h = fct_display(uint16(IMAGE),round(2.54/Resolution));
            elseif strcmp(ans,'90 left')
                tIMAGE(:,:,1) = flipud(IMAGE(:,:,1)');
                tIMAGE(:,:,2) = flipud(IMAGE(:,:,2)');
                tIMAGE(:,:,3) = flipud(IMAGE(:,:,3)');
                IMAGE = tIMAGE; clear tIMAGE;
                close(h);figure('NumberTitle','off','Name','Averaged image');
                h = fct_display(uint16(IMAGE),round(2.54/Resolution));
            end
            close(h);figure('NumberTitle','off','Name','Averaged image');
            h = fct_display(uint16(IMAGE),round(2.54/Resolution));
            ans = questdlg('Do you want to flip the image horizontally?','Image','Left-right','Up-down','No','No') ;
            if strcmp(ans,'Left-right')
                tIMAGE(:,:,1) = fliplr(IMAGE(:,:,1));
                tIMAGE(:,:,2) = fliplr(IMAGE(:,:,2));
                tIMAGE(:,:,3) = fliplr(IMAGE(:,:,3));
                IMAGE = tIMAGE; clear tIMAGE;
                close(h);figure('NumberTitle','off','Name','Averaged image');
                h = fct_display(uint16(IMAGE),round(2.54/Resolution));
            elseif strcmp(ans,'Up-down')
                tIMAGE(:,:,1) = flipud(IMAGE(:,:,1));
                tIMAGE(:,:,2) = flipud(IMAGE(:,:,2));
                tIMAGE(:,:,3) = flipud(IMAGE(:,:,3));
                IMAGE = tIMAGE; clear tIMAGE;
                close(h);figure('NumberTitle','off','Name','Averaged image');
                h = fct_display(uint16(IMAGE),round(2.54/Resolution));
            end
            ans = questdlg('Do you want to modify the image again?','Image','Yes','No','No') ;
            if strcmp(ans,'No')
                flg=0;
            end
            close(h);
        end
        flg = 1;
        %%%%%
        ans = 'No';%skip the following
        %ans = questdlg('Do you want to calibrate the signal?','Image','Yes','No','No') ;
        if strcmp(ans,'Yes')
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
            IMAGE(:,:,1) = uint16((2^16-1)*interp1(Sr,I,double(IMAGE(:,:,1))/(2^16-1),'pchip'));
            IMAGE(:,:,2) = uint16((2^16-1)*interp1(Sg,I,double(IMAGE(:,:,2))/(2^16-1),'pchip'));
            IMAGE(:,:,3) = uint16((2^16-1)*interp1(Sb,I,double(IMAGE(:,:,3))/(2^16-1),'pchip'));
            IMAGE = uint16(IMAGE);
            tfilename = tfilename(:);
            tfilename = [tfilename(1:end-7)' '_ave_lin.tif'];
        else
            tfilename = tfilename(:);
            tfilename = [tfilename(1:end-7)' '_ave.tif'];           
        end
        [ofilename,opathname] = uiputfile({'*.tif'},'Save average film',tfilename);
        if ofilename==0
        else
            ofilename = fct_makecleanfilename(opathname,ofilename);
            imwrite(uint16(IMAGE),ofilename,'tif','Resolution',round(2.54/Resolution));
        end
    end
end

% --------------------------------------------------------------------
function GRIDTRIGGER_Callback(hObject, eventdata, handles)
% hObject    handle to gridTRIGGER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.gridflag == 0
    handles.gridflag = 1;
else
    handles.gridflag = 0;
    handles.gridres = fct_autoruler(handles.z,handles.DELTA);
end

%UPDATE DISPLAY
handles = fct_updatedisplay(handles);;

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function DIST_TOOL_Callback(hObject, eventdata, handles)
% hObject    handle to DIST_TOOL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    imdistline;
end

% --------------------------------------------------------------------
function PROF_VERTICAL_Callback(hObject, eventdata, handles)
% hObject    handle to PROF_VERTICAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    button = questdlg('What selection type?','ROI','Free','Rectangular','Point','Free');
    [x,profil]=fct_getmeanprofilefromselect(handles.z,handles.DELTA,'y',button,1,[]);
    if length(x)~=0
        fct_showprofile(x,profil,handles.DELTA);
        button = questdlg('Do you want to save profile in ascii format?','Profile','Yes','No','No');
        if strcmp(button,'Yes')
            [ifilename,ipathname]=uiputfile({'*.plot'},'Save profile');
            ACTUALDIR=cd;
            filename = fct_makecleanfilename(ipathname,ifilename);
            xtmp = x(:) - (max(x)+min(x))/2;
            ytmp = profil(:);
            A = cat(2,xtmp,ytmp);
            fid = fopen(filename,'w');
            for ii=1:size(A,1)
                fprintf(fid,'%f %f\n',A(ii,1),A(ii,2));
            end
            fclose(fid);
        end        
    end
end

% --------------------------------------------------------------------
function PROF_HORIZONTAL_Callback(hObject, eventdata, handles)
% hObject    handle to PROF_HORIZONTAL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    button = questdlg('What selection type?','ROI','Free','Rectangular','Point','Free');
    [x,profil] = fct_getmeanprofilefromselect(handles.z,handles.DELTA,'x',button,1,[]);
    if length(x)~=0
       fct_showprofile(x,profil,handles.DELTA);
        button = questdlg('Do you want to save profile in ascii format?','Profile','Yes','No','No');
        if strcmp(button,'Yes')
            [ifilename,ipathname]=uiputfile({'*.plot'},'Save profile');
            ACTUALDIR=cd;
            filename = fct_makecleanfilename(ipathname,ifilename);
            xtmp = x(:) - (max(x)+min(x))/2;
            ytmp = profil(:);
            A = cat(2,xtmp,ytmp);
            fid = fopen(filename,'w');
            for ii=1:size(A,1)
                fprintf(fid,'%f %f\n',A(ii,1),A(ii,2));
            end
            fclose(fid);
        end
    end
end
% --------------------------------------------------------------------
function MAP_Callback(hObject, eventdata, handles)
% hObject    handle to RULER5mm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    answer = 'Yes';
    if handles.filter==0
        answer = questdlg('Data was not filtered. Continue?','Filtering','Yes','No','No');
    end
    if strcmp(answer,'Yes')
        button = questdlg('How do you wish to normalize?','Normalization','Point','Region','Maximum value','Region');
        if strcmp(button,'Point')
            [x,y,nx,ny] = fct_getpoints(handles.z,handles.DELTA);
            z = handles.z;
            norm = z(ny,nx);
        elseif strcmp(button,'Region')
            figure(handles.H);
            [A,center,owidth] = fct_getroi(handles.z,handles.DELTA,'Rectangular',0);
            norm = mean(mean(A));
        elseif strcmp(button,'Maximum value')
            norm = max(max(handles.z));
        end
        
        [x,y,nx,ny]=fct_getpoints(handles.z,handles.DELTA);
        
        z = handles.z*100/norm;
        ox=mean(x);
        oy=mean(y);
        xx=((1:size(z,2))-(size(z,2)+1)/2)*handles.DELTA-ox;
        yy=((1:size(z,1))-(size(z,1)+1)/2)*handles.DELTA-oy;
        [X,Y]=meshgrid(xx,yy);
        
        %DISPLAY
        figure;
        maxi = max(max(z));
        mini = floor(min(min(z)));
        %N=[10 20 30 40 50 60 70 80 85 90 95 100 105 110 115 maxi];
        increment = max(1,floor((maxi-mini)/10));
        N = mini:increment:maxi;
        [C,h]=contourf(X,Y,z,N);
        clabel(C,h,'Fontsize',14);
        view(0,-90);
        set(h,'Linewidth',2);
        set(gca,'Linewidth',1);
        %set(gca,'ytick',0:0.2:10,'xtick',-10:0.5:10);
        set(gca,'XLimMode','manual');
        set(gca,'YLimMode','manual');
        set(gca,'Xgrid','on');
        set(gca,'Ygrid','on');
        axis image;
        set(gca,'XLim',[min(xx) max(xx)]);
        set(gca,'YLim',[min(yy) max(yy)]);
        colormap(handles.color);
        %shading('interp');
    end
    
end

% --------------------------------------------------------------------
function RPVERT_Callback(hObject, eventdata, handles)
% hObject    handle to PROFVERT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    button = questdlg('What selection type?','ROI','Free','Rectangular','Point','Free');
    [y,rp] = fct_getmeanprofilefromselect(handles.z,handles.DELTA,'y',button,1,[]);
    if length(y)~=0
        y = y -y(1);
        fct_showrp(y,rp,handles.DELTA);
        button = questdlg('Do you want to save profile in ascii format?','Profile','Yes','No','No');
        if strcmp(button,'Yes')
            [ifilename,ipathname]=uiputfile({'*.plot'},'Save profile');
            ACTUALDIR=cd;
            filename = fct_makecleanfilename(ipathname,ifilename);
            xtmp = y(:) - (max(y)+min(y))/2;
            ytmp = rp(:);
            A = cat(2,xtmp,ytmp);
            fid = fopen(filename,'w');
            for ii=1:size(A,1)
                fprintf(fid,'%f %f\n',A(ii,1),A(ii,2));
            end
            fclose(fid);
        end
    end
end

% --------------------------------------------------------------------
function RPHORI_Callback(hObject, eventdata, handles)
% hObject    handle to MEANTIF16RGB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    button = questdlg('What selection type?','ROI','Free','Rectangular','Point','Free');
    [x,rp] = fct_getmeanprofilefromselect(handles.z,handles.DELTA,'x',button,1,[]);
    if length(x)~=0
        x = x-x(1);
        fct_showrp(x,rp,handles.DELTA);        
        button = questdlg('Do you want to save profile in ascii format?','Profile','Yes','No','No');
        if strcmp(button,'Yes')
            [ifilename,ipathname]=uiputfile({'*.plot'},'Save profile');
            ACTUALDIR=cd;
            filename = fct_makecleanfilename(ipathname,ifilename);
            xtmp = x(:) - (max(x)+min(x))/2;
            ytmp = rp(:);
            A = cat(2,xtmp,ytmp);
            fid = fopen(filename,'w');
            for ii=1:size(A,1)
                fprintf(fid,'%f %f\n',A(ii,1),A(ii,2));
            end
            fclose(fid);
        end
    end
end

% --------------------------------------------------------------------
function ISODOSES_Callback(hObject, eventdata, handles)
% hObject    handle to gridTRIGGER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    answer = 'Yes';
    if handles.filter==0
        answer = questdlg('Data was not filtered. Continue?','Filtering','Yes','No','No');
    end
    if strcmp(answer,'Yes')
        button = questdlg('How do you wish to normalize?','Normalizsation','Point','Region','Maximum value','Region');
        if strcmp(button,'Point')
            [x,y,nx,ny] = fct_getpoints(handles.z,handles.DELTA);
            z = handles.z;
            norm = z(ny,nx);
        elseif strcmp(button,'Region')
            figure(handles.H);
            [A,center,owidth] = fct_getroi(handles.z,handles.DELTA,'Rectangular',0);
            norm = mean(mean(A));
        elseif strcmp(button,'Maximum value')
            norm = max(max(handles.z));
        end
        
        [x,y,nx,ny]= fct_getpoints(handles.z,handles.DELTA);
        
        z = handles.z/norm*100;
        ox=mean(x);
        oy=mean(y);
        xx=((1:size(z,2))-(size(z,2)+1)/2)*handles.DELTA-ox;
        yy=((1:size(z,1))-(size(z,1)+1)/2)*handles.DELTA-oy;
        [X,Y]=meshgrid(xx,yy);
        
        %DISPLAY
        figure;
        maxi = max(max(z));
        mini = max(1,floor(min(min(z))));
        %N=[10 20 30 40 50 60 70 80 85 90 95 100 105 110 115 maxi];
        increment = floor((maxi-mini)/10);
        N = mini:increment:maxi;
        [C,h]=contour(X,Y,z,N);
        clabel(C,h,'Fontsize',14);
        view(0,-90);
        set(h,'Linewidth',3);
        set(gca,'Linewidth',1);
        %set(gca,'ytick',0:0.2:10,'xtick',-10:0.5:10);
        set(gca,'XLimMode','manual');
        set(gca,'YLimMode','manual');
        set(gca,'Xgrid','on');
        set(gca,'Ygrid','on');
        axis image;
        set(gca,'XLim',[min(xx) max(xx)]);
        set(gca,'YLim',[min(yy) max(yy)]);
        colormap(handles.color);
    end
    
end

% --------------------------------------------------------------------
function EXPORT_TIF16RGB_Callback(hObject, eventdata, handles)
% hObject    handle to EXPORT_TIF16RGB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    z = handles.z;
    h = handles.H;
    flag = 0;
    if (handles.data==1)
        flag = 1;
    elseif (handles.data==2)
        z = 10.^(-z).*(2^(handles.BITSini)-1);
        flag = 1;
    elseif (handles.data==3)
        msgbox('Cannot export .tif file from net optical density.','Error');
    elseif (handles.data==4)
        msgbox('Cannot export .tif file from dose.','Error');
    elseif (handles.data==5)
        z = 10.^(-z).*(2^(handles.BITSini)-1);
        flag = 1;
    end
    if flag
        [ofilename,opathname] = uiputfile({'*.tif'},'.tif image export');
        if ofilename==0            
        else
            filename = fct_makecleanfilename(opathname,ofilename);
            A = uint16(z);
            imwrite(A(:,:,1),filename,'tif','Resolution',round(2.54/handles.DELTA));
        end
    else
        h_msg = gcf;
        figure(handles.H);
        handles = fct_updatedisplay(handles);;
        figure(h_msg);
    end
    
end

% --------------------------------------------------------------------
function HOMOG_CREATECORRMAT_Callback(hObject, eventdata, handles)
% hObject    handle to HOMOG_CREATECORRMAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles = fct_CreateHomogCorrectionMatrixAuto(handles);
% handles = fct_CreateHomogCorrectionMatrixSemiAuto(handles);
if 1%CCD-based
    handles = fct_CreateHomogCorrectionMatrixDetailed(handles);
else%automated/continuous
    handles = fct_CreateHomogCorrectionMatrixAutomated(handles);
end
% --------------------------------------------------------------------
function HOMOG_CORRECT_Callback(hObject, eventdata, handles)
% hObject    handle to HOMOG_CORRECT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[ifilename1,ipathname1] = uigetfile({'*.tif'},'Choose .tif image to be corrected','multiselect','on');

handles = fct_updatedisplay(handles);
figure(handles.H);

if iscell(ifilename1)
    nbfilms = numel(ifilename1);
else
    nbfilms =1;
    tmp = ifilename1;
    clear ifilename1;
    ifilename1{1} = tmp;
end
if tmp==0
    return;
end
[cfilename,cpathname] = uigetfile({'*.hm3'},'Choose correction matrix file');        
if (~strcmp(class(cfilename),'double'))
    flag = 1;
    corrfname = fct_makecleanfilename(cpathname,cfilename);
    [rawres,rawimsize,corrdir,xkrange,rrange,grange,brange,pr,pg,pb] = fct_ReadCorrmatrixDetailed(corrfname);
else
    flag = 0;
end
if flag
    for ii=1:nbfilms
        if (~strcmp(class(ifilename1),'double'))

            [Im,res,bits,channel] = fct_read_tif_image(fct_makecleanfilename(ipathname1,ifilename1{ii}),'All');
            [x,y] = fct_gridindextopos(size(Im,1),size(Im,2),res);
            
            if channel==5
                figure;
                h = fct_display(uint16(Im),res);

%                 direction = questdlg('In which direction do you correct the image?','Direction','Horizontal','Vertical','Horizontal') ;
%                 switch direction
%                     case 'Horizontal'
%                         corrdir = 2;
%                         x = x;
%                     case 'Vertical'
%                         corrdir = 1;
%                         x = y;
%                 end
%                 [Imcorr,xcorr] = fct_CorrectHomogAutomated(x,Im,PR,PG,PB,xunique,corrdir);

                [Imcorr,xkcorr] = fct_CorrectHomogDetailed(Im,rawres,rawimsize,corrdir,xkrange,rrange,grange,brange,pr,pg,pb);
                
                close(h);

                figure;
                h = fct_display(uint16(Imcorr),res);

                tmp = ifilename1{ii};
                defaultname = sprintf('%s_homog.tif',tmp(1:(length(tmp)-4)));
                [ofilename,opathname] = uiputfile({'*.tif'},'Save corrected image',defaultname);

                if (~strcmp(class(ofilename),'double'))
                    filename = fct_makecleanfilename(opathname,ofilename);
                    imwrite(Imcorr,filename,'tif','Resolution',round(2.54/res));
                end
                
                close(h);
            else
                errordlg('Homogeneity correction applies to 16 RGB images only.')
            end
        end
    end
end
% --------------------------------------------------------------------
function HOMOG_VISUALIZE_Callback(hObject, eventdata, handles)
% hObject    handle to HOMOG_VISUALIZE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ifilename,ipathname] = uigetfile({'*.hm3'},'Choose correction matrix file');

%handles = fct_updatedisplay(handles);;
%figure(handles.H);

if (~strcmp(class(ifilename),'double'))
    
    fname = fct_makecleanfilename(ipathname,ifilename);
    if 1%detailed
        [corrres,rawimsize,corrdir,xkrange,rrange,grange,brange,pr,pg,pb] = fct_ReadCorrmatrixDetailed(fname);
        [xgrid,ygrid] = fct_gridindextopos(rawimsize(1),rawimsize(2),corrres);
        if corrdir==1
            xrange = ygrid(xkrange(1):xkrange(2));
        elseif ocrrdir==2
            xrange = xgrid(xkrange(1):xkrange(2));
        end
        % Vizualise correction map   
        sr = [min(rrange):100:max(rrange)]';
        sg = [min(grange):100:max(grange)]';
        sb = [min(brange):100:max(brange)]';
        Sr = repmat(sr,1,length(xrange));
        Sg = repmat(sg,1,length(xrange));
        Sb = repmat(sb,1,length(xrange));
        Xr = repmat(xrange(:)',length(sr),1);
        Xg = repmat(xrange(:)',length(sg),1);
        Xb = repmat(xrange(:)',length(sb),1);
        PR1 = repmat(pr(1,:),length(sr),1);
        PR2 = repmat(pr(2,:),length(sr),1);
        PR3 = repmat(pr(3,:),length(sr),1);
        PG1 = repmat(pg(1,:),length(sg),1);
        PG2 = repmat(pg(2,:),length(sg),1);
        PG3 = repmat(pg(3,:),length(sg),1); 
        PB1 = repmat(pb(1,:),length(sb),1);
        PB2 = repmat(pb(2,:),length(sb),1);
        PB3 = repmat(pb(3,:),length(sb),1);  
        Rhat = min(65535,max(0,PR1.*Sr.^2+PR2.*Sr.^1+PR3.*Sr.^0));
        Ghat = min(65535,max(0,PG1.*Sg.^2+PG2.*Sg.^1+PG3.*Sg.^0));
        Bhat = min(65535,max(0,PB1.*Sb.^2+PB2.*Sb.^1+PB3.*Sb.^0));

        figure;
        subplot(1,3,1);
        mesh(Xr,Sr,Rhat./Sr);
        title('Red');
        subplot(1,3,2);
        mesh(Xg,Sg,Ghat./Sg);
        title('Green');
        subplot(1,3,3);
        mesh(Xb,Sb,Bhat./Sb);
        title('Blue');
    else%automated/continuous
        [PR,PG,PB,xunique,rlims,glims,blims]  = fct_ReadHomogCorrAutomated(fname);
        %here we don't care about limits in signal because we take the signal
        %correction as linear
        %%%
        [IX,IS] = meshgrid(xunique,0:100:65535);
        [lin,col] = size(IX); 
        ix = reshape(IX,lin*col,1);
        is = reshape(IS,lin*col,1);
        F = @(x,s,n) cat(2,repmat(s(:),1,n+1).^repmat(zeros(1,n+1),length(s(:)),1).*repmat(x(:),1,n+1).^repmat(0:n,length(x(:)),1),...
                       repmat(s(:),1,n+1).^repmat(ones(1,n+1) ,length(s(:)),1).*repmat(x(:),1,n+1).^repmat(0:n,length(x(:)),1)    );

        orderR = length(PR)/2-1;
        orderG = length(PG)/2-1;
        orderB = length(PB)/2-1;
        %predicted correction factors for interpolation

        [IXr,ISr] = meshgrid(xunique,rlims(1):100:rlims(2));
        [IXg,ISg] = meshgrid(xunique,glims(1):100:glims(2));
        [IXb,ISb] = meshgrid(xunique,blims(1):100:blims(2));
        [linr,colr] = size(IXr); 
        [ling,colg] = size(IXg); 
        [linb,colb] = size(IXb); 
        ixr = reshape(IXr,linr*colr,1);
        ixg = reshape(IXg,ling*colg,1);
        ixb = reshape(IXb,linb*colb,1);
        isr = reshape(ISr,linr*colr,1);
        isg = reshape(ISg,ling*colg,1);
        isb = reshape(ISb,linb*colb,1);
        %predicted correction factors for interpolation
        Fr = reshape(F(ixr,isr,orderR)*PR,linr,colr)./ISr;
        Fg = reshape(F(ixg,isg,orderG)*PG,ling,colg)./ISg;
        Fb = reshape(F(ixb,isb,orderB)*PB,linb,colb)./ISb;
        %Display of raw signal inference
        figure;
        subplot(1,3,1);
        mesh(IXr,ISr,Fr); colormap('jet')
        view(-15,15);
        title('Red');
        subplot(1,3,2);
        mesh(IXg,ISg,Fg); colormap('jet')
        view(-15,15); 
        title('Green');
        subplot(1,3,3);
        mesh(IXb,ISb,Fb); colormap('jet')
        view(-15,15);
        title('Blue');
    end      
    
end

% --------------------------------------------------------------------
function HOMOG_MATRXDIFF_Callback(hObject, eventdata, handles)
% hObject    handle to HOMOG_MATRXDIFF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This is now diabled

if 0
    [ifilename,ipathname] = uigetfile({'*.hm3'},'Choose correction matrix file');

    %handles = fct_updatedisplay(handles);;
    %figure(handles.H);

    if (~strcmp(class(ifilename),'double'))

        fname = fct_makecleanfilename(ipathname,ifilename);
    %     [xx,s1,s2,s3,c1,c2,c3,p1,p2,p3,order1,order2,sym,m,n,bits] = fct_ReadCorrmatrixAuto(fname);

        [xi,y1i,y2i,y3i,c1i,c2i,c3i,x,y1,y2,y3,c1,c2,c3] = fct_ReadCorrmatrixSemiAuto(fname);

        c1hat = interp2(xi,y1i,c1i,x,y1,'spline');
        c2hat = interp2(xi,y2i,c2i,x,y2,'spline');
        c3hat = interp2(xi,y3i,c3i,x,y3,'spline');    
    %     F1 = []; F2 = []; F3 = [];
    %     for i = 0:sym:order1
    %         m1 = []; m2 = []; m3 = [];
    %         for j=0:1:order2
    %             m1 = cat(2,m1,xx.^i.*s1.^j);
    %             m2 = cat(2,m2,xx.^i.*s2.^j);
    %             m3 = cat(2,m3,xx.^i.*s3.^j);
    %         end
    %         F1 = cat(2,F1,m1);
    %         F2 = cat(2,F2,m2);
    %         F3 = cat(2,F3,m3);
    %     end
    % 
    %     c1hat = reshape(F1*p1,m,n);
    %     c2hat = reshape(F2*p2,m,n);
    %     c3hat = reshape(F3*p3,m,n);
    % 
    %     s = sprintf('Max order: in x = %d, in OD = %d',order1,order2);
    % 
    %     x = reshape(xx,m,n);
    %     y1 = reshape(s1,m,n);y2 = reshape(s2,m,n);y3 = reshape(s3,m,n);
    %     c1 = reshape(c1,m,n); c2 = reshape(c2,m,n); c3 = reshape(c3,m,n);
    %     
    %     s = sprintf('Max orders (position,signal) = (%d,%d)',order1,order2);
    % 

        h1 = figure; 
    %     subplot(1,3,1);mesh(x,y1,c1./c1hat-1);xlabel('Position');ylabel('Signal');zlabel('Relative error');
    %     title('Red');
    %     subplot(1,3,2);mesh(x,y2,c2./c2hat-1);xlabel('Position');ylabel('Signal');zlabel('Relative error');
    %     title('Green');
    %     subplot(1,3,3);mesh(x,y3,c3./c3hat-1);xlabel('Position');ylabel('Signal');zlabel('Relative error');
    %     title('Blue');
        subplot(1,3,1);plot3(x(:),y1(:),c1(:)./c1hat(:)-1,'.r','markersize',1);xlabel('Position');ylabel('Signal');zlabel('Relative error');
        title('Red'); grid('on');
        subplot(1,3,2);plot3(x(:),y2(:),c2(:)./c2hat(:)-1,'.g','markersize',1);xlabel('Position');ylabel('Signal');zlabel('Relative error');
        title('Green'); grid('on');
        subplot(1,3,3);plot3(x(:),y3(:),c3(:)./c3hat(:)-1,'.b','markersize',1);xlabel('Position');ylabel('Signal');zlabel('Relative error');
        title('Blue'); grid('on');
    %     subplot(2,3,4);plot(reshape(x,m*n,1),reshape(c1-c1hat,m*n,1),'.r');xlabel('Position'); ylabel('Red residual');
    %     subplot(2,3,5);plot(reshape(x,m*n,1),reshape(c2-c2hat,m*n,1),'.g');xlabel('Position'); ylabel('Green residual');
    %     subplot(2,3,6);plot(reshape(x,m*n,1),reshape(c3-c3hat,m*n,1),'.b');xlabel('Position'); ylabel('Blue residual');
    end
    % [ifilename1,ipathname1] = uigetfile({'*.hmg'},'Choose correction matrix  #1');
    % button = questdlg('Which method do you wish to use?','Method','Polynomial','Filter','New','Filter') ;
    % if strcmp(button,'Polynomial')
    %     method1 = 1;
    % elseif strcmp(button,'Filter')
    %     method1 = 0;
    % else
    %     method1 = 2;
    % end
    % [ifilename2,ipathname2] = uigetfile({'*.hmg'},'Choose correction matrix  #2');
    % button = questdlg('Which method do you wish to use?','Method','Polynomial','Filter','New','Filter') ;
    % if strcmp(button,'Polynomial')
    %     method2 = 1;
    % elseif strcmp(button,'Filter')
    %     method2 = 0;
    % else
    %     method2 = 2;
    % end
    % 
    % handles = fct_updatedisplay(handles);;
    % figure(handles.H);
    % 
    % if (~strcmp(class(ifilename1),'double'))&&(~strcmp(class(ifilename2),'double'))
    %     
    %     fname1 = fct_makecleanfilename(ipathname1,ifilename1);
    %     if (method1==1) %method polynoms
    %         [p,center,xmin,xmax,smin,smax] = fct_readcorrmatrix(fname1);
    %         [P,XI,SI] = fct_makecorrmatrix(p,center,[xmin xmax],[smin smax]);
    %         xlim1 = [xmin xmax];
    %         slim1 = [smin smax];
    %     elseif (method1==0)%method filter
    %         [N,xi,F,center] = fct_readcorrmatrix2(fname1);
    %         [P,XI,SI] = fct_makecorrmatrix2(N,xi,F,[min(xi) max(xi)],[min(center) max(center)]);
    %         xlim1 = [min(xi) max(xi)];
    %         slim1 = [min(center) max(center)];
    %     elseif (method1==2)%method new
    %         [p,xmin,xmax,x0,smin,smax] = fct_readcorrmatrix3(fname1);
    %         [P,XI,SI] = fct_makecorrmatrix3(p,[xmin xmax],x0,[smin smax]);
    %         xlim1 = [xmin xmax];
    %         slim1 = [smin smax];
    %     end
    %     
    %     X1 = XI;
    %     S1 = SI;
    %     P1 = P;
    %     
    %     fname2 = fct_makecleanfilename(ipathname2,ifilename2);
    %     if (method2==1) %method polynoms
    %         [p,center,xmin,xmax,smin,smax] = fct_readcorrmatrix(fname2);
    %         [P,XI,SI] = fct_makecorrmatrix(p,center,[xmin xmax],[smin smax]);
    %         xlim2 = [xmin xmax];
    %         slim2 = [smin smax];
    %     elseif (method2==0) %method filter
    %         [N,xi,F,center] = fct_readcorrmatrix2(fname2);
    %         [P,XI,SI] = fct_makecorrmatrix2(N,xi,F,[min(xi) max(xi)],[min(center) max(center)]);
    %         xlim2 = [min(xi) max(xi)];
    %         slim2 = [min(center) max(center)];
    %     elseif (method2==2) %method new
    %         [p,xmin,xmax,x0,smin,smax] = fct_readcorrmatrix3(fname2);
    %         [P,XI,SI] = fct_makecorrmatrix3(p,[xmin xmax],x0,[smin smax]);
    %         xlim2 = [xmin xmax];
    %         slim2 = [smin smax];
    %     end
    %     
    %     X2 = XI;
    %     S2 = SI;
    %     P2 = P;
    %     
    %     xlim = [max(xlim1(1),xlim2(1)) min(xlim1(2),xlim2(2))];
    %     slim = [max(slim1(1),slim2(1)) min(slim1(2),slim2(2))];
    %     
    %     Nbin = 128;
    %     xmin = min(xlim);
    %     xmax = max(xlim);
    %     deltax = (xmax-xmin)/Nbin;
    %     xi = xmin:deltax:xmax;
    %     smin = min(slim);
    %     smax = max(slim);
    %     deltas = (smax-smin)/Nbin;
    %     si = smin:deltas:smax;
    %     
    %     [XI,SI] = meshgrid(xi,si);
    %     
    %     P1I = interp2(X1,S1,P1,XI,SI,'linear');
    %     P2I = interp2(X2,S2,P2,XI,SI,'linear');
    %     h2 = figure('NumberTitle','off','Name','Difference of correction matrices');
    %     mesh(SI,XI,P2I-P1I);
    %     
    %     ylabel('Off-axis diatance (cm)');
    %     xlabel('Signal');
    %     zlabel('Factor difference');
    %     
    %     D = (P2I-P1I)./(P1I+P2I)*2;
    %     %assumes relative uncertainty is constant
    %     rms = sqrt(sum(sum(D.*D))/size(D,1)/size(D,2));
    %     D = log(65535./(P2I.*SI))/log(10) - log(65535./(P1I.*SI))/log(10);
    %     bias = sum(sum(D))/size(D,1)/size(D,2);
    %     %Lors que l'on compare une methode filtree (1 cm^2) avec la m�thode polynomes, cette difference
    %     %devrait etre de l'ordre de l'incertitude obtenue avec la courbe de calibration.
    %     dmb = sprintf('OD rms difference = %.4f\nOD bias = %.4f',rms/log(10),bias);
    %     msgbox(dmb,'Statistical analysis');
    % end
end
% --------------------------------------------------------------------
function GROUP_PIX_Callback(hObject, eventdata, handles)
% hObject    handle to GROUP_PIX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)    
    N  = inputdlg({'Entrer the factor by which you want to reduce the resolution in each dimension:'},'Pixel grouping',1,{'1'}) ;
    if max(size(N))==0
    else
        N = str2double(N);
        N = max(1,round(N));
        handles.z = fct_reducematrix(handles.z,N);
        handles.DELTA = handles.DELTA*N;
        handles = fct_updatedisplay(handles);;
        % Update handles structure
        guidata(hObject, handles);
    end
end


% --------------------------------------------------------------------
function DEFINE_CHANNEL_Callback(hObject, eventdata, handles)
% hObject    handle to DEFINE_CHANNEL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~fct_isthereanimage(handles)
    
%     button = questdlg('Which channel do you want to analyze?','Channel','Red','Green','Blue','Red');
    str{1} = 'Red';
    str{2} = 'Green';
    str{3} = 'Blue';
    str{4} = 'Multi';
    [i,ok] = listdlg('Name','Channel','ListString',str);
    if (ok==1)&&(length(i)==1)
        if length(i)==0
            i = 1;
        end
        handles.defaultcolor =  str{i};
    end
    
    handles = fct_updatedisplay(handles);;
    % Update handles structure
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function REPEATED_MEASUREMENTS_Callback(hObject, eventdata, handles)
% hObject    handle to REPEATED_MEASUREMENTS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.channel==4
    
else
%     button = questdlg('How do you wish to enter normalization values?','Values','File','Manual','File') ;
%     DOSE = [];
%     if strcmp(button,'Manual')
%         nb  = inputdlg({'Number of values'},'Number of values',1);
%         nb  = str2double(nb);
%         nprompt = ceil(nb/10);
%         for j = 1:nprompt
%             clear prompt;
%             for i= 1+10*(j-1):min(10*j,nb)
%                 dumb = sprintf('Film #%d\nDose (CMU) or NOD',i);
%                 prompt{i-10*(j-1)} = dumb;
%             end
%             answer = inputdlg(prompt,'Dose values (CMU) or NOD',1);
%             DOSE = cat(1,DOSE,str2double(answer));
%         end
%         DOSE = DOSE';
%         %ici il devrait y avoir une facon de valider les valeurs de dose
%         DOSE = abs(DOSE);
%     else
%         [ifilename,ipathname] = uigetfile({'*.txt'},'Choose file containing normalization dose or NOD values');
%         if ~strcmp(class(ifilename),'double')
%             file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');
%             DOSE = fscanf(file,'%f',[1 inf]);
%         else
%             for j = 1:nprompt
%                 clear prompt;
%                 for i= 1+10*(j-1):min(10*j,nb)
%                     dumb = sprintf('Film #%d\nDose (CMU) or NOD',i);
%                     prompt{i-10*(j-1)} = dumb;
%                 end
%                 answer = inputdlg(prompt,'Dose values (CMU) or NOD',1);
%                 DOSE = cat(1,DOSE,str2double(answer));
%             end
%             DOSE = DOSE';
%         end
%     end
%     DOSE1 = DOSE;
%     MU1  = inputdlg({'Number of MU'},'Number of MU',1);
%     MU1  = str2double(MU1);
% 
%     button = questdlg('How do you wish to enter values to be normalized?','Values','File','Manual','File') ;
%     DOSE = [];
%     if strcmp(button,'Manual')
%         nb  = inputdlg({'Number of values'},'Number of values',1);
%         nb  = str2double(nb);
%         nprompt = ceil(nb/10);
%         for j = 1:nprompt
%             clear prompt;
%             for i= 1+10*(j-1):min(10*j,nb)
%                 dumb = sprintf('Film #%d\nDose (CMU) or NOD',i);
%                 prompt{i-10*(j-1)} = dumb;
%             end
%             answer = inputdlg(prompt,'Dose values (CMU) or NOD',1);
%             DOSE = cat(1,DOSE,str2double(answer));
%         end
%         DOSE = DOSE';
%         %ici il devrait y avoir une facon de valider les valeurs de dose
%         DOSE = abs(DOSE);
%     else
%         [ifilename,ipathname] = uigetfile({'*.txt'},'Choose file containing dose or NOD values');
%         if ~strcmp(class(ifilename),'double')
%             file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');
%             DOSE = fscanf(file,'%f',[1 inf]);
%         else
%             for j = 1:nprompt
%                 clear prompt;
%                 for i= 1+10*(j-1):min(10*j,nb)
%                     dumb = sprintf('Film #%d\nDose (CMU) or NOD',i);
%                     prompt{i-10*(j-1)} = dumb;
%                 end
%                 answer = inputdlg(prompt,'Dose values (CMU) or NOD',1);
%                 DOSE = cat(1,DOSE,str2double(answer));
%             end
%             DOSE = DOSE';
%         end
%     end
%     DOSE2 = DOSE;
%     MU2  = inputdlg({'Number of MU'},'Number of MU',1);
%     MU2  = str2double(MU2);
% 
%     [ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');
%     filename = [ipathname ifilename]
%     file = fopen(filename,'r');
%     [DOSE,OD,M,type,sigparam,Npix] = fct_readcalfile(file);
%     fclose(file);
%     [odi,dosei] = fct_getdose_from_od(DOSE,OD,M,type);
% 
%     if max(DOSE1)<1 %net OD
%         DOSE1 = interp1(odi,dosei,DOSE1,'linear');
%     end
%     if max(DOSE2)<1 %net OD
%         DOSE2 = interp1(odi,dosei,DOSE2,'linear');
%     end
% 
%     N1 =length(DOSE1);
%     N2 =length(DOSE2);
% 
%     dose1 = mean(DOSE1);
%     V = fct_getcovarmatrix(DOSE1(:),Npix(1),Npix(1),filename);
%     sdose1 = dose1*sqrt( sum(sum(V))/N1^2/dose1^2);
% 
%     [dose1 sdose1]
% 
%     dose2 = mean(DOSE2);
%     V = fct_getcovarmatrix(DOSE2(:),Npix(1),Npix(1),filename);
%     sdose2 = dose2*sqrt( sum(sum(V))/N2^2/dose2^2);
% 
%     [dose2 sdose2]
% 
%     dose_vect = [DOSE1(:)' DOSE2(:)']';
%     V = fct_getcovarmatrix(dose_vect,Npix(1),Npix(1),filename);
%     sd0 = sum(sum(V(1:N1,1:N1)))/(dose1^2*N1*N2);
%     sd1 = sum(sum(V((N1+1):(N1+N2),(N1+1):(N1+N2))))/(dose1^2*N1*N2);
%     sd2 = sum(sum(V(1:N1,(N1+1):(N1+N2))))/(dose1*dose2*N1*N2);
%     sd = sqrt(sd0+sd1-2*sd2 + 2*0.001^2/(N1*N2));
%     rdf = dose2/dose1*MU1/MU2;
%     srdf = sd*rdf;
% 
%     [rdf srdf]
% 
%     dumb=sprintf('Dose #1 = %.1f +/- %.1f\nDose #2 = %.1f +/- %.1f\nRatio = %.4f +/- %.4f\n\n*This does not include type B uncertainties.',dose1,sdose1,dose2,sdose2,rdf,srdf);
%     msgbox(dumb,'Results');
%     h = gcf;
%     figure(handles.H);
%     handles = fct_updatedisplay(handles);;
%     figure(h);
end


% --------------------------------------------------------------------
function STACKIMAGES_Callback(hObject, eventdata, handles)
% hObject    handle to STACKIMAGES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ifilename1,ipathname1] = uigetfile({'*.tif'},'Choose .tif file #1 to be stacked');
if ifilename1~=0
    [ifilename2,ipathname2] = uigetfile({'*.tif'},'Choose .tif file #2 to be stacked');
    if ifilename2~=0
        if (~strcmp(class(ifilename1),'double'))&&(~strcmp(class(ifilename2),'doube'))
            cpu = computer;
            if strcmp(cpu,'PCWIN')||strcmp(cpu,'PCWIN64')
                c = '\'; % dos
            else
                c = '/'; % mac/linux/unix
                %In the case someone uses a SUN, I don't know if this is correct
            end
            path = ipathname1;
            filename = ifilename1;
            if (path(max(size(path,1),size(path,2)))==c)
                name = [path filename];
            else
                name = [path c filename];
            end
            ifilename1 = name;

            path = ipathname2;
            filename = ifilename2;
            if (path(max(size(path,1),size(path,2)))==c)
                name = [path filename];
            else
                name = [path c filename];
            end
            ifilename2 = name;

            [IMAGE1,Resolution1,bits1,channel1] = fct_read_tif_image(ifilename1,'All');
            [IMAGE2,Resolution2,bits2,channel2] = fct_read_tif_image(ifilename2,'All');
            s1 = size(IMAGE1);
            s2 = size(IMAGE2);
            if ((Resolution1~=Resolution2)||(bits1~=bits2)||(channel1~=channel2))
                errordlg('Images are not compatible.');
            elseif (s1(1)~=s2(1))||(s1(2)~=s2(2))
                errordlg('Images are not compatible.');
            else
                Resolution = Resolution1;
                channel = channel1;
                button = questdlg('In which direction you wish to stack the images?','Direction','Vertical','Horizontal','Horizontal');
                if(strcmp(button,'Vertical'))
                    if channel==5
                        x = [IMAGE1(:,:,1)' IMAGE2(:,:,1)']'; %[a b] = cat(2,a,b) while [a' b']' = cat(1,a,b)
                        y = [IMAGE1(:,:,2)' IMAGE2(:,:,2)']';
                        z = [IMAGE1(:,:,3)' IMAGE2(:,:,3)']';
                    elseif channel==4
                        x = [IMAGE1' IMAGE2']';
                    end
                else
                    if channel==5
                        x = [IMAGE1(:,:,1) IMAGE2(:,:,1)];
                        y = [IMAGE1(:,:,2) IMAGE2(:,:,2)];
                        z = [IMAGE1(:,:,3) IMAGE2(:,:,3)];
                    else channel==4
                        x = [IMAGE1 IMAGE2];
                    end 
                end
                [N,M] = size(x);
                if channel==5
                    NEWIMAGE = uint16(zeros(N,M,3));
                    NEWIMAGE(:,:,1) = x;
                    NEWIMAGE(:,:,2) = y;
                    NEWIMAGE(:,:,3) = z;        
                elseif channel==4
                    NEWIMAGE = uint16(x);
                end
                [ofilename,opathname,FilterIndex] = uiputfile({'*.tif'},'Save stacked image');
                if ofilename==0
                else
                    path = opathname;
                    filename = ofilename;
                    if (path(max(size(path,1),size(path,2)))==c)
                        name = [path filename];
                    else
                        name = [path c filename];
                    end
                    ofilename = name; 
                    imwrite(NEWIMAGE,ofilename,'tif','Resolution',round(2.54/Resolution));        
                end
            end
        end
    end
end
% --------------------------------------------------------------------
function MERGESPLITIMAGES_Callback(hObject, eventdata, handles)
% hObject    handle to MERGESPLITIMAGES (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ifilename1,ipathname1] = uigetfile({'*.tif'},'Choose .tif file #1 to be merged');
if ifilename1~=0
    [ifilename2,ipathname2] = uigetfile({'*.tif'},'Choose .tif file #2 to be merged');
    if ifilename2~=0

        if (~strcmp(class(ifilename1),'double'))&&(~strcmp(class(ifilename2),'double'))

            cpu = computer;
            if strcmp(cpu,'PCWIN')||strcmp(cpu,'PCWIN64')
                c = '\'; % dos
            else
                c = '/'; % mac/linux/unix
                %In the case someone uses a SUN, I don't know if this is correct
            end

            path = ipathname1;
            filename = ifilename1;
            if (path(max(size(path,1),size(path,2)))==c)
                name = [path filename];
            else
                name = [path c filename];
            end
            ifilename1 = name;

            path = ipathname2;
            filename = ifilename2;
            if (path(max(size(path,1),size(path,2)))==c)
                name = [path filename];
            else
                name = [path c filename];
            end
            ifilename2 = name;

            [IMAGE1,Resolution1,bits1,channel1] = fct_read_tif_image(ifilename1,'All');
            [IMAGE2,Resolution2,bits2,channel2] = fct_read_tif_image(ifilename2,'All');
            s1 = size(IMAGE1);
            s2 = size(IMAGE2);
            if ((Resolution1~=Resolution2)||(bits1~=bits2)||(channel1~=channel2))
                errordlg('Images are not compatible.');
            elseif (s1(1)~=s2(1))||(s1(2)~=s2(2))
                errordlg('Images are not compatible.');
            else
                Resolution = Resolution1;
                channel = channel1;
                button = questdlg('In which direction you wish to merge the images?','Direction','Vertical','Horizontal','Horizontal');

                if(strcmp(button,'Vertical'))
                    index = s1(1)/2;
                    if (index-floor(index))==0
                        jump = 1;
                        vx = [];
                        vy = [];
                        vz = [];
                    else
                        jump = 2;
                        index = index - 0.5;
                        if channel==5
                            vx = IMAGE1(index+1,:,1)/2 + IMAGE1(index+1,:,1)/2;
                            vy = IMAGE1(index+1,:,2)/2 + IMAGE1(index+1,:,2)/2;
                            vz = IMAGE1(index+1,:,3)/2 + IMAGE1(index+1,:,3)/2;
                        elseif channel==4
                            vx = IMAGE1(index+1,:)/2 + IMAGE1(index+1,:)/2;
                        end
                    end
                    if channel==5
                        x = [IMAGE1(1:index,:,1)' vx' IMAGE2(index+jump:s1(1),:,1)']'; %[a b] = cat(2,a,b) while [a' b']' = cat(1,a,b)
                        y = [IMAGE1(1:index,:,2)' vy' IMAGE2(index+jump:s1(1),:,2)']';
                        z = [IMAGE1(1:index,:,3)' vz' IMAGE2(index+jump:s1(1),:,3)']';
                    elseif channel==4
                        x = [IMAGE1(1:index,:)' vx' IMAGE2(index+jump:s1(1),:)']'; %[a b] = cat(2,a,b) while [a' b']' = cat(1,a,b)
                    end
                else
                    index = s1(2)/2;
                    if (index-floor(index))==0
                        jump = 1;
                        vx = [];
                        vy = [];
                        vz = [];
                    else
                        jump = 2;
                        index = index - 0.5;
                        if channel==5
                            vx = IMAGE1(:,index+1,1)/2 + IMAGE1(:,index+1,1)/2;
                            vy = IMAGE1(:,index+1,2)/2 + IMAGE1(:,index+1,2)/2;
                            vz = IMAGE1(:,index+1,3)/2 + IMAGE1(:,index+1,3)/2;
                        elseif channel==4
                            vx = IMAGE1(:,index+1)/2 + IMAGE1(:,index+1)/2;    
                        end
                    end
                    if channel==5
                        x = [IMAGE1(:,1:index,1) vx IMAGE2(:,index+jump:s1(2),1)];
                        y = [IMAGE1(:,1:index,2) vy IMAGE2(:,index+jump:s1(2),2)];
                        z = [IMAGE1(:,1:index,3) vz IMAGE2(:,index+jump:s1(2),3)];
                    elseif channel==4
                        x = [IMAGE1(:,1:index) vx IMAGE2(:,index+jump:s1(2))];
                    end
                end

                [N,M] = size(x);
                if channel==5
                    NEWIMAGE = uint16(zeros(N,M,3));
                    NEWIMAGE(:,:,1) = x;
                    NEWIMAGE(:,:,2) = y;
                    NEWIMAGE(:,:,3) = z;
                elseif channel==4
                    NEWIMAGE = uint16(x); 
                end
                [ofilename,opathname,FilterIndex] = uiputfile({'*.tif'},'Save merged image');
                path = opathname;
                filename = ofilename;
                if (path(max(size(path,1),size(path,2)))==c)
                    name = [path filename];
                else
                    name = [path c filename];
                end
                ofilename = name;

                imwrite(NEWIMAGE,ofilename,'tif','Resolution',round(2.54/Resolution));

            end
        end
    end
end
% --------------------------------------------------------------------
function VARANALYSIS_SINGLE_Callback(hObject, eventdata, handles)
% hObject    handle to VARANALYSIS_SINGLE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if strcmp(handles.defaultcolor,'Multi')
%     [ifilename,ipathname] = uigetfile({'*.mlt'},'Choose calibration curve');
% else
%     [ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');
% end
[ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');

if ~strcmp(class(ifilename),'double')
    
    filename = fct_makecleanfilename(ipathname,ifilename);
    if 0%strcmp(handles.defaultcolor,'Multi')
        file = fopen(filename,'r');
        [DOSE,XI,sXI,c,V,res,sximesh] = fct_ReadCalFileMulti(file);
        fclose(file);       
        npix = floor(0.1/handles.CCDres^2);
        dose = 10:10:max(DOSE);
        G = [dose(:).^0 dose(:).^1];
        xi = G*c;
        sd1 = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,sximesh.npix0,dose);
        sd2 = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,npix,dose);

        k = find((sd2(:)./dose(:))==min(sd2./dose(:)));
        r = 0.1:0.01:2;
        sr1 = fct_UncertaintyMultichannelSingleMeasRel(c,V,sximesh,sximesh.npix0,dose(k),sximesh.npix0,r(:).*dose(k));
        sr2 = fct_UncertaintyMultichannelSingleMeasRel(c,V,sximesh,npix,dose(k),npix,r(:).*dose(k));
        
        figure;
        subplot(1,2,1)
        plot(dose(:),sd1(:)./dose(:)*100,'--k',dose(:),sd2(:)./dose(:)*100,'k','linewidth',2);
        set(gca,'Fontsize',12,'Fontweight','bold','Fontname','Times New Roman');
        xlabel('DI');
        ylabel('Type A uncertainty on DI (%)');
        title(fct_addbackslash(ifilename));
        set(gca,'Ylim',[0 5],'Xlim',[0 max(DOSE)]);
        grid on;
        str{1} = sprintf('ROI = %.2f X %.2f mm^2',res*10,res*10);
        str{2} = sprintf('ROI = 1 mm^2');
        legend(str,'Location','northeast');

        [dose(:) sd2(:)./dose(:)*100]
        
        subplot(1,2,2);
%         plot(r(:),sr1(:)./r(:)*100,'k',r(:),sr2(:)./r(:)*100,'r',r(:),sr(:)./r(:)*100,'b','linewidth',2);
        plot(r(:),sr1(:)./r(:)*100,'--k',r(:),sr2(:)./r(:)*100,'k','linewidth',2);
        set(gca,'Fontsize',12,'Fontweight','bold','Fontname','Times New Roman');
        xlabel('DI ratio');
        ylabel('Type A uncertainy on DI ratio (%)');
        title(fct_addbackslash(ifilename));
        set(gca,'Ylim',[0 5],'Xlim',[0 max(r)]);
        grid on;
        str{1} = sprintf('ROI = %.2f X %.2f mm^2',res*10,res*10);
        str{2} = sprintf('ROI = 1 mm^2');
        legend(str,'Location','northeast');
        
        [r(:) sr2(:)./r(:)*100]
    else
        file = fopen(filename,'r');
        [DOSE,OD,M,type,sigparam,Npix] = fct_readcalfile(file);
        fclose(file);

        dose= 10:10:max(DOSE);
        V = fct_getcovarmatrix(dose,Npix(1),Npix(2),filename);
        sd = sqrt(diag(V));

        figure;
        plot(dose(:),sd(:)./dose(:)*100,'k','linewidth',2);
        set(gca,'Fontsize',12,'Fontweight','bold','Fontname','Times New Roman');
        xlabel('Dose (CMU)');
        ylabel('Relative uncertainty (%)');
        title(fct_addbackslash(ifilename));
        set(gca,'Ylim',[0 5],'Xlim',[0 max(DOSE)]);
        grid on;

        [dose(:) sd(:)./dose(:)*100]

        k= find((sd(:)./dose(:))==min(sd(:)./dose(:)));

        dosenorm = dose(k)
        r = dose(:)/dosenorm;
        sr = r*0;
        k = length(dose);
        for i=1:k
            V = fct_getcovarmatrix([dose(i) dosenorm],Npix(1),Npix(2),filename);
            sr(i)    = r(i)*sqrt(V(1,1)/dose(i)^2+V(2,2)/dosenorm^2 - 2*V(1,2)/dose(i)/dosenorm);
        end

        [r(:) sr(:)./r(:)*100]

        figure;
        plot(r(:),sr(:)./r(:)*100,'k','linewidth',2);
        set(gca,'Fontsize',12,'Fontweight','bold','Fontname','Times New Roman');
        xlabel('Dose ratio');
        ylabel('Relative uncertainty (%)');
        title(fct_addbackslash(ifilename));
        set(gca,'Ylim',[0 5],'Xlim',[0 max(r)]);
        grid on;
    end    
    
end

% --------------------------------------------------------------------
function VARANALYSIS_REPEATED_Callback(hObject, eventdata, handles)
% hObject    handle to VARANALYSIS_REPEATED (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if strcmp(handles.defaultcolor,'Multi')
%     [ifilename,ipathname] = uigetfile({'*.mlt'},'Choose calibration curve');
% else
%     [ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');
% end

[ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');

if ~strcmp(class(ifilename),'double')
    
    filename = fct_makecleanfilename(ipathname,ifilename);
    if 0%strcmp(handles.defaultcolor,'Multi')
        file = fopen(filename,'r');
        [DOSE,XI,sXI,c,V,res,sximesh] = fct_ReadCalFileMulti(file);
        fclose(file);

        npix = floor(0.1/handles.CCDres^2);
        dose = 10:10:max(DOSE);
        G = [dose(:).^0 dose(:).^1];
        xi = G*c;
        sd1 = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,sximesh.npix0,dose);
        sd2 = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,npix,dose);

        k = find((sd2(:)./dose(:))==min(sd2./dose(:)));
                       
        N = 10;
        sdose = zeros(N,1);
        sr = zeros(N,1);
        for i=1:N
            sdose(i) = fct_UncertaintyMultichannelMultipleMeas(c,V,sximesh,npix,dose(k),i);
            sr(i) = fct_UncertaintyMultichannelMultipleMeasRel(c,V,res,sximesh,dose(k),npix,dose(k),npix,i);
        end
        
        figure;
        subplot(1,2,1);
        semilogx((1:N)',sdose(:)./dose(k)*100,'k','linewidth',2);
        set(gca,'Fontsize',12,'Fontweight','bold','Fontname','Times New Roman');
        xlabel('Number of measurements');
        ylabel('Type A uncertainty on DI (%)');
        title(fct_addbackslash(ifilename));
        set(gca,'Ylim',[0 max(sdose(:)./dose(k))*100*1.1],'Xlim',[1 N]);
        grid on;

        subplot(1,2,2);
        semilogx((1:N)',sr*100,'k','linewidth',2);
        set(gca,'Fontsize',12,'Fontweight','bold','Fontname','Times New Roman');
        xlabel('Number of measurements');
        ylabel('Type A uncertainty on DI ratio (%)');
        title(fct_addbackslash(ifilename));
        set(gca,'Ylim',[0 max(sr)*100*1.1],'Xlim',[1 N]);
        grid on;

    else
        file = fopen(filename,'r');
        [DOSE,OD,M,type,sigparam,Npix] = fct_readcalfile(file);
        fclose(file);

        dose= 1:1:max(DOSE);
        V = fct_getcovarmatrix(dose,Npix(1),Npix(2),filename);
        sd = sqrt(diag(V));
        k= find((sd(:)./dose(:))==min(sd(:)./dose(:)));

        dosemin = dose(k)

        N = 100;
        DOSE1 = ones(N,1)*dosemin;

        SDOSE = zeros(N,1);

        for i=1:N
            V = fct_getcovarmatrix(DOSE1(1:i),Npix(1),Npix(2),filename);
            SDOSE(i) = sqrt( sum(sum(V))/i^2/dosemin^2);
        end

        [(1:N)' SDOSE(:)*100]

        figure;
        semilogx((1:N)',SDOSE*100,'k','linewidth',2);
        set(gca,'Fontsize',12,'Fontweight','bold','Fontname','Times New Roman');
        xlabel('Number of measurements');
        ylabel('Absolute dose relative uncertainty (%)');
        title(fct_addbackslash(ifilename));
        set(gca,'Ylim',[0 max(SDOSE)*100*1.1],'Xlim',[1 N]);
        grid on;

        DOSE2 = DOSE1;

        SR = zeros(N,1);

        dose_vect = [DOSE1(:)' DOSE2(:)']';
        for i=1:N
            V = fct_getcovarmatrix(dose_vect,Npix(1),Npix(1),filename);
            sd0 = sum(sum(V(1:i,1:i)))/(dosemin^2*i*i);
            sd1 = sum(sum(V((i+1):(i+i),(i+1):(i+i))))/(dosemin^2*i*i);
            sd2 = sum(sum(V(1:i,(i+1):(i+i))))/(dosemin*dosemin*i*i);
            SR(i) = sqrt(sd0+sd1-2*sd2 + 2*0.001^2/(i*i));
        end

        [(1:N)' SR(:)*100]

        figure;
        semilogx((1:N)',SR*100,'k','linewidth',2);
        set(gca,'Fontsize',12,'Fontweight','bold','Fontname','Times New Roman');
        xlabel('Number of measurements');
        ylabel('Dose ratio relative uncertainty (%)');
        title(fct_addbackslash(ifilename));
        set(gca,'Ylim',[0 max(SR)*100*1.1],'Xlim',[1 N]);
        grid on;
    end
    
end


% --------------------------------------------------------------------
function VARANALYSIS_UNCERT_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to VARANALYSIS_UNCERT_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if strcmp(handles.defaultcolor,'Multi')
%     [ifilename,ipathname] = uigetfile({'*.mlt'},'Choose calibration curve');
% else
%     [ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');
% end

[ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');

if ~strcmp(class(ifilename),'double')
    
    filename = fct_makecleanfilename(ipathname,ifilename);
    file = fopen(filename,'r');
        if 0%strcmp(handles.defaultcolor,'Multi')
            file = fopen(filename,'r');
            [DOSE,XI,sXI,c,V,rescal,sximesh] = fct_ReadCalFileMulti(file);
            fclose(file);

            figure;
            mesh(sximesh.x,sximesh.y,sximesh.z);
            xlabel('sqrt(N_0/N_1)');
            ylabel('xi');
            zlabel('sxi');

%             dose = 1:1:max(DOSE);
%             sd = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,sximesh.npix0,dose);
%             k = find((sd(:)./dose(:))==min(sd./dose(:)));
%                        
%             ROI = (0.1:0.1:100)/10; %in cm
%             n = (ROI/handles.CCDres).^2;
%             G = [dose(k).^0 dose(k).^1];
%             xi = G*c;
%             sd = ROI*0;
%             sxi = ROI*0;
%             for i=1:length(sd)
%                 sd(i) = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,n(i),dose(k));
%                 sxi(i) = fct_GetsXI(sximesh,n(i),xi);
%             end
% 
%             xisymb = char(958);
%             figure; 
%             subplot(1,2,1);
%             semilogx(ROI*10,sxi,'k','Linewidth',2);
%             ylabel([ 'Type A uncertainty on ' xisymb],'Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
%             xlabel('ROI width (mm)','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
%             set(gca,'Xlim',[min(ROI*10) max(ROI*10)],'Ylim',[0 max(sxi)*1.1]);
%             grid on;
%             title(fct_addbackslash((ifilename)));    
%             
%             subplot(1,2,2);
%             semilogx(ROI*10,sd./dose(k)*100,'k','Linewidth',2);
%             ylabel('Type A uncertainty on DI (%)','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
%             xlabel('ROI width (mm)','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
%             set(gca,'Xlim',[min(ROI*10) max(ROI*10)],'Ylim',[0 max(sd./dose(k)*100)*1.1]);
%             grid on;
%             title(fct_addbackslash((ifilename))); 
%             
        else
            [DOSE,OD,M,opt,sigparam,Npix] = fct_readcalfile(file);
            fclose(file);
            [p,sy,R,df] = fct_lsf(DOSE,OD,M,opt);
            ROI = (0.1:0.1:100)/10; %in cm
            N = (ROI/handles.CCDres).^2;
            ROI = ROI*10; %in mm
            s = sqrt(2*sigparam(1).^2+(1./N+1/Npix(2))*sigparam(2).^2);

            figure;
            semilogx(ROI,s,'k','Linewidth',2);
            ylabel('NOD uncertainty','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
            xlabel('ROI (mm)','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
            set(gca,'Xlim',[min(ROI) max(ROI)],'Ylim',[0 max(s)*1.1]);
            grid on;
            title(fct_addbackslash((ifilename)));
        end
end


% --------------------------------------------------------------------
function VARANALYSIS_RESIDUALS_Callback(hObject, eventdata, handles)
% hObject    handle to VARANALYSIS_RESIDUALS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if strcmp(handles.defaultcolor,'Multi')
%     [ifilename,ipathname] = uigetfile({'*.mlt'},'Choose calibration curve');
% else
%     [ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');
% end

[ifilename,ipathname] = uigetfile({'*.cal'},'Choose calibration curve');

if ~strcmp(class(ifilename),'double')

    filename = fct_makecleanfilename(ipathname,ifilename);
    file = fopen(filename,'r');
     if 0%strcmp(handles.defaultcolor,'Multi')
        file = fopen(filename,'r');
        [DOSE,XI,sXI,c,V,rescal,sximesh] = fct_ReadCalFileMulti(file);
        fclose(file);

        npix = floor(1/handles.CCDres^2);
        G = [DOSE(:).^0 DOSE(:).^1];
        xi = G*c;
        sxihat = fct_GetsXI(sximesh,npix,XI);
        dose = (XI -c(1))/c(2);
        sdose = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,npix,dose);
        
        DOSE = DOSE(:)
        dose = dose(:)
%         [DOSE (dose-DOSE)./DOSE*100 sdose./DOSE*100]
        k = find(DOSE~=0);

        xisymb = char(958);
        figure;
        subplot(1,2,1);
        errorbar(DOSE(k), (dose(k)-DOSE(k))./DOSE(k)*100,2* sdose(k)./DOSE(k)*100,'ko','Linewidth',2);
        ylabel('Error (%)','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
        xlabel('DI','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
%         set(gca,'Xlim',[0 max(DOSE)+min(DOSE)]);
        grid on;
        title(fct_addbackslash(ifilename));
        legend('Calibration (k=2)','Location','northeast');
        
        subplot(1,2,2);
        errorbar(xi, xi-XI, 2*sxihat,'ko','Linewidth',2);
        ylabel('Error','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
        xlabel(xisymb,'Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
        set(gca,'Xlim',[-0.1*max(xi) max(xi)*1.1]);
        grid on;
        title(fct_addbackslash(ifilename));
        legend('Calibration (k=2)','Location','northwest');
     else
        [DOSE,OD,M,opt,sigparam,Npix] = fct_readcalfile(file);
        fclose(file);
        [odi,dosei] = fct_getcalcurvepoints(DOSE,OD,M,opt);
        dose = interp1(odi,dosei,OD);
        k = intersect(find(~isnan(dose)),find(DOSE~=0));
        dose = dose(k);
        DOSE = DOSE(k);
        V = fct_getcovarmatrix(dose,Npix(1),Npix(2),filename);

        DOSE = DOSE(:);
        dose = dose(:);
        sdose = sqrt(diag(V));
        sdose = sdose(:);
%         [DOSE (dose-DOSE)./DOSE*100 sdose./DOSE*100]

        figure;
        errorbar(DOSE, (dose-DOSE)./DOSE*100, sdose./DOSE*100,'ko','Linewidth',2);
        ylabel('Error (%)','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
        xlabel('Dose (CMU)','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
        set(gca,'Xlim',[0 max(DOSE)+min(DOSE)]);
        grid on;
        title(fct_addbackslash(ifilename));

        file = fopen(filename,'r');
        [DOSE,OD,M,opt,sigparam,Npix] = fct_readcalfile(file);
        fclose(file);
        [p,sy,R,df] = fct_lsf(DOSE,OD,M,opt);
        F = fct_F_matrix(DOSE,M,opt);
        ODfit = F*p;
        r = OD(:)-ODfit(:);
        sr = sy+0*r;

        figure;
        errorbar(ODfit, r, sr,'ko','Linewidth',2);
        ylabel('Error','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
        xlabel('Net optical density','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
        set(gca,'Xlim',[-0.1*max(ODfit) max(ODfit)*1.1]);
        grid on;
        title(fct_addbackslash(ifilename));
     end
end

% --------------------------------------------------------------------
function CLOSEALL_Callback(hObject, eventdata, handles)
% hObject    handle to CLOSEALL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if fct_isthereanimage(handles)
    
    clear handles.H;

    handles = fct_initGafgui(handles.defaultcolor);
%     %Initialization of main variables
%     handles.dataini = 0;
%     handles.SATMINini = 0;
%     handles.SATMAXini = 0;
%     handles.BCKGRNDini = [];
%     handles.NpixBCKGRNDini = 0;
%     handles.DOSEREFini = [];
%     handles.NpixDOSEREFini = [];
%     handles.DELTAini = 0;
%     handles.channelini = fct_colortochannel(handles.defaultcolor);
%     handles.filterini = 0;
%     handles.ORIGINini = [0 0 ]';
%     handles.Z = [];
%     handles.cfilenameini = '';
%     handles.BITSini = 0;
%     handles.data = handles.dataini;
%     handles.SATMIN = handles.SATMINini;
%     handles.SATMAX = handles.SATMAXini;
%     handles.BCKGRND = handles.BCKGRNDini;
%     handles.NpixBCKGRND = handles.NpixBCKGRNDini;
%     handles.DOSEREF = handles.DOSEREFini;
%     handles.NpixDOSEREF = handles.NpixDOSEREFini;
%     handles.DELTA = handles.DELTAini;
%     handles.channel = handles.channelini;
%     handles.filter = handles.filterini;
%     handles.ORIGIN = handles.ORIGINini;
%     handles.z = handles.Z;
%     handles.cfilename = handles.cfilenameini;
    
    %Initialization of other variables
    if handles.channel==4
        handles.color = 'gray';
    else
        handles.color ='jet';
    end
    handles.gridflag = 0;
    handles.gridres = 1;
    handles.ifilename = '';
    handles = fct_updatedisplay(handles);
    
    % Update handles structure
    guidata(hObject, handles);
end

% --------------------------------------------------------------------
function MENULINEAR_Callback(hObject, eventdata, handles)
% hObject    handle to MENULINEAR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function LINEAR_CORRECT_Callback(hObject, eventdata, handles)
% hObject    handle to LINEAR_CORRECT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[ifilename1,ipathname1] = uigetfile({'*.tif'},'Choose .tif image to be corrected','multiselect','on');

handles = fct_updatedisplay(handles);
figure(handles.H);

if iscell(ifilename1)
    nbfilms = numel(ifilename1);
else
    nbfilms =1;
    tmp = ifilename1;
    clear ifilename1;
    ifilename1{1} = tmp;
end
if tmp==0
    return;
end
[cfilename,cpathname] = uigetfile({'*.lin'},'.lin signal calibration');
if (~strcmp(class(cfilename),'double'))
    flag = 1;
    filename = fct_makecleanfilename(cpathname,cfilename);
    LUT  = fct_ReadSignalCorrection(filename)
else
    flag = 0;
end
if flag
    for ii=1:nbfilms
        if (~strcmp(class(ifilename1),'double'))

            [Im,res,bits,channel] = fct_read_tif_image(fct_makecleanfilename(ipathname1,ifilename1{ii}),'All');

            if channel==5
                
                R = double(Im(:,:,1));
                G = double(Im(:,:,2));
                B = double(Im(:,:,3));
                [lin,col,~] = size(Im);
                
%                 R = reshape(interp1(double(LUT.raw),double(LUT.red),  reshape(R,lin*col,1)),lin,col);
%                 G = reshape(interp1(double(LUT.raw),double(LUT.green),reshape(G,lin*col,1)),lin,col);
%                 B = reshape(interp1(double(LUT.raw),double(LUT.blue), reshape(G,lin*col,1)),lin,col);
                
                R = 65535*fct_SignaltoTransmission(Im(:,:,1),LUT,1);
                G = 65535*fct_SignaltoTransmission(Im(:,:,2),LUT,2);
                B = 65535*fct_SignaltoTransmission(Im(:,:,3),LUT,3);
                
                clear Imcorr;
                
                Imcorr(:,:,1) = uint16(R);
                Imcorr(:,:,2) = uint16(G);
                Imcorr(:,:,3) = uint16(B);
                
                tmp = ifilename1{ii};
                defaultname = sprintf('%s_lin.tif',tmp(1:(length(tmp)-4)));
                [ofilename,opathname] = uiputfile({'*.tif'},'Save corrected image',defaultname);

                if (~strcmp(class(ofilename),'double'))
                    filename = fct_makecleanfilename(opathname,ofilename);
                    imwrite(Imcorr,filename,'tif','Resolution',round(2.54/res));
                end
                
            else
                errordlg('Homogeneity correction applies to 16 RGB images only.')
            end
        end
    end
end

% --------------------------------------------------------------------
function LINEAR_CREATE_Callback(hObject, eventdata, handles)
% hObject    handle to LINEAR_CREATE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = fct_CreateSignalLinearizationMap(handles);

% --------------------------------------------------------------------
function LINEAR_VISUALIZE_Callback(hObject, eventdata, handles)
% hObject    handle to LINEAR_VISUALIZE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ofilename,opathname] = uigetfile({'*.lin'},'.lin signal calibration');
if ofilename==0            
else
    filename = fct_makecleanfilename(opathname,ofilename);
    LUT = fct_ReadSignalCorrection(filename);
    figure;
    hold on;
    plot(LUT.raw,LUT.red,'-r','linewidth',2); 
    plot(LUT.raw,LUT.green,'-g','linewidth',2); 
    plot(LUT.raw,LUT.blue,'-b','linewidth',2);
    legend('Red channel - LUT','Green channel - LUT','Blue channel - LUT','location','northwest');
    hold off;
    xlabel('Raw signal');
    ylabel('Corrected signal');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------------------------------------------------------------------------------------------------------------------------------%
% FUNCTIONS
% ----------------------------------------------------------------------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Functions are (now) outside Gafgui in Gafgui_functions folder

