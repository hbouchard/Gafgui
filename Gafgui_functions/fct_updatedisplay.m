% --------------------------------------------------------------------
function handles = fct_updatedisplay(handles)

%this chooses the background color
color = fct_ModeToColorVect(fct_colortochannel(handles.defaultcolor));

%this updates the menu so that only the valid ones are enabled
handles = fct_ShowMenus(handles);

if min(size(handles.z))==0
    himage = imagesc(handles.z);
    set(gcf,'Color',color);
    axis('off');
    %pixval (gca,'off');
else
    nlines = size(handles.z,1);
    ncols = size(handles.z,2);
    [x,y] = fct_gridindextopos(nlines,ncols,handles.DELTA);
    x = x - handles.ORIGIN(1);
    y = y - handles.ORIGIN(2);
    axis('on');
    himage = imagesc(x,y,handles.z);
    set(gcf,'Color',color);
    vect1 = [handles.gridres*ceil(min(x)/handles.gridres)-handles.gridres:handles.gridres:handles.gridres*floor(max(x)/handles.gridres)+handles.gridres];
    vect2 = [handles.gridres*ceil(min(y)/handles.gridres)-handles.gridres:handles.gridres:handles.gridres*floor(max(y)/handles.gridres)+handles.gridres];
    set(gca,'DataAspectRatio',[1 1 1],'Xtick',vect1,'Ytick',vect2);
    set(gcf,'Units','pixels');
    %flipud gray map
    if strcmp(handles.color,'gray')
        if handles.data==1
            colormap(handles.color);
        else
            c = colormap(handles.color);
            colormap(flipud(c));
        end
    else
        colormap(handles.color);
    end
    htext = impixelinfoval(gcf,himage);
    %set(htext,'FontWeight','bold');
    set(htext,'FontSize',10);
    %pixval(gca,'on');
    colorbar('vert');
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %Jasmine Duchaine 30 March 2020: fix bug with _ 
    % set(gca,'defaulttextinterpreter','latex'); 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    title((fct_guiinfo(handles)),'interpreter','none');
    if handles.gridflag == 1
        grid on;
    else
        grid off;
    end
end
handles.H = gcf;
