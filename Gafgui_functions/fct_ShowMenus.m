function handles = fct_ShowMenus(handles);
%%%%
handles
if isfield(handles, 'menus')
  if isfield(handles.menus, 'MENUFILE')
%this confirms that we are in Gafgui...
    %Menu1 (File): some should remain on at all times
    set(handles.menus.MENUFILE, 'Enable', 'on');
        set(handles.menus.MENUIMPORTIMAGE, 'Enable', 'on');
            set(handles.menus.IMPORT_MATLAB, 'Enable', 'on');
            set(handles.menus.IMPORT_TIF, 'Enable', 'on');
        set(handles.menus.MENU_SAVE_IMAGE, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
            set(handles.menus.SAVE_MATLAB, 'Enable', 'on');
        set(handles.menus.MENUEXPORT, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
            set(handles.menus.EXPORT_TIF16RGB, 'Enable', 'on');
        set(handles.menus.CLOSEALL, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
        set(handles.menus.QUIT, 'Enable', 'on');
    %Menu2 (Tools): those should remain on only when there is an image
    set(handles.menus.MENUTOOLS, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
        set(handles.menus.RESTOREIM, 'Enable', 'on');
        set(handles.menus.MENUAFF, 'Enable', 'on');
            set(handles.menus.GRIDTRIGGER, 'Enable', 'on');
            set(handles.menus.GRIDSIZE, 'Enable', 'on');
            set(handles.menus.COLORMAP, 'Enable', 'on');
        set(handles.menus.ROISELECT, 'Enable', 'on');
        set(handles.menus.MENUROT, 'Enable', 'on');
            set(handles.menus.ROTRIGHT, 'Enable', 'on');        
            set(handles.menus.ROTLEFT, 'Enable', 'on');
            set(handles.menus.MIRVERT, 'Enable', 'on');
            set(handles.menus.MIRHORIZ, 'Enable', 'on');
            set(handles.menus.ROTPOINTS, 'Enable', 'on');
        set(handles.menus.DIST_TOOL, 'Enable', 'on');
        set(handles.menus.GROUP_PIX, 'Enable', 'on');
        set(handles.menus.FOLD_IMAGE, 'Enable', 'on');
        set(handles.menus.MENUFILTERIM, 'Enable', fct_LogicToSwitch(handles.filter==0));
            set(handles.menus.FILTER_RECT, 'Enable', 'on');
            set(handles.menus.FILTER_GAUSS, 'Enable', 'on');    
        set(handles.menus.MENU_SATURATE, 'Enable', 'on');
        set(handles.menus.SATURATE_MIN, 'Enable', 'on');
        set(handles.menus.SATURATE_MAX, 'Enable', 'on');
    %Menu3 (Preparation): those should remain on when there is no image
    set(handles.menus.MENUPREP, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)~=1));
        set(handles.menus.MEANTIF16RGB, 'Enable', 'on');
        set(handles.menus.STACKIMAGES, 'Enable', 'on');
        set(handles.menus.MERGESPLITIMAGES, 'Enable', 'on');
        set(handles.menus.HOMOG_CORRECT, 'Enable', 'on');
        set(handles.menus.MULTICHANNEL_CORRECT, 'Enable', 'on');
    %Menu4 (Data): those depend on the data type
    set(handles.menus.MENUDONNE, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
        set(handles.menus.SHOWSIGNAL, 'Enable', fct_LogicToSwitch((handles.data==1)||(handles.data==2)));
        set(handles.menus.SHOWRAWOD, 'Enable', fct_LogicToSwitch((handles.data==1)||(handles.data==2)));
        set(handles.menus.BCKGRND_DEFINE, 'Enable', fct_LogicToSwitch(handles.data==2));
        set(handles.menus.SHOWNETOD, 'Enable', fct_LogicToSwitch(((handles.data==2)||(handles.data==3))&&(numel(handles.BCKGRND)~=0)));
 %       set(handles.menus.SHOWXI, 'Enable', fct_LogicToSwitch(handles.data==5));
        set(handles.menus.SHOWDOSE, 'Enable', fct_LogicToSwitch((handles.data==3)||(handles.data==4)));
%         set(handles.menus.DEFINE_REF_DOSE, 'Enable', fct_LogicToSwitch(handles.data==4));
        set(handles.menus.DEFINE_REF_DOSE, 'Enable', 'off');%I am switching this off temporarily 20 jan 2021 
    %Menu5 (Analysis): those should remain on with image ony
    set(handles.menus.MENUANALYZE, 'Enable', 'on');
        set(handles.menus.MENU_MEANDOSE, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
            set(handles.menus.MEANDOSE_SELECT, 'Enable', 'on');
            set(handles.menus.MEANDOSE_MULTIPLE, 'Enable', 'on');
        set(handles.menus.MEAN_MEANDOSE_ORIGIN, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
            set(handles.menus.MEANDOSE_ORIGIN, 'Enable', 'on');
            set(handles.menus.SET_ORIGIN, 'Enable', 'on');
        set(handles.menus.MENUREGISTRATION, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
            set(handles.menus.REGISTRATION_CROSSHAIR, 'Enable', 'on');
%         set(handles.menus.REPEATED_MEASUREMENTS, 'Enable', fct_LogicToSwitch(handles.channel==4));
        set(handles.menus.REPEATED_MEASUREMENTS, 'Enable', 'off'); %I am switching this off temporarily 17 jan 2021 
        set(handles.menus.MENUPROFILE, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
            set(handles.menus.PROF_VERTICAL, 'Enable', 'on');
            set(handles.menus.PROF_HORIZONTAL, 'Enable', 'on');
        set(handles.menus.MENURP, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
            set(handles.menus.RPVERT, 'Enable', 'on');
            set(handles.menus.RPHORI, 'Enable', 'on');
        set(handles.menus.MENU2DDIST, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
            set(handles.menus.ISODOSES, 'Enable', 'on');
            set(handles.menus.MAP, 'Enable', 'on');
    %Menu6 (Characterization): this should remain on at all times
    set(handles.menus.MENUCHARACT, 'Enable', 'on');
        set(handles.menus.MENUHOMOG, 'Enable', 'on');
            set(handles.menus.HOMOG_CREATECORRMAT, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)~=1));
            set(handles.menus.HOMOG_VISUALIZE, 'Enable', 'on');
            set(handles.menus.HOMOG_MATRXDIFF, 'Enable', 'off');%I am switching this off temporarily 17 jan 2020 (1 year priori to the paper)
        set(handles.menus.MENULINEAR, 'Enable', 'on');
            set(handles.menus.LINEAR_CREATE,'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)~=1));
        set(handles.menus.MULTICHANNEL_CHARACTERIZE, 'Enable', 'on');
            set(handles.menus.MULTI_CREATE, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)~=1));
            set(handles.menus.MULTI_VISUALIZE, 'Enable', 'off'); %I am switching this off for now 20 jan 2021
        set(handles.menus.MENUCALIBRATION, 'Enable', 'on');
            set(handles.menus.CALCURVE_SETDEFAULT, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)==1));
            set(handles.menus.CALCURVE_CREATE, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)~=1));
            set(handles.menus.CALCURVE_VISUALIZE, 'Enable', 'on');    
    %Menu7 (Variance analysis): this should remain on at all times
%     set(handles.menus.MENU_VARANALYSIS, 'Enable', 'on');
    set(handles.menus.MENU_VARANALYSIS, 'Enable', 'off');%I am switching this off temporarily 20 jan 2021 
        set(handles.menus.VARANALYSIS_SINGLE, 'Enable', 'on');
        set(handles.menus.VARANALYSIS_REPEATED, 'Enable', 'on');
        set(handles.menus.VARANALYSIS_UNCERT_ROI, 'Enable', 'on');
        set(handles.menus.VARANALYSIS_RESIDUALS, 'Enable', 'on');
    %Menu8 (Parameters): this should remain on if no image only
    set(handles.menus.MENUPARAM, 'Enable', fct_LogicToSwitch(fct_isthereanimage(handles)~=1));
        set(handles.menus.DEFINE_CHANNEL, 'Enable', 'on');
    %Menu9 (?): this should remain on at all times
    set(handles.menus.MENUABOUT, 'Enable', 'on');
        set(handles.menus.VERSION, 'Enable', 'on');
  end
end