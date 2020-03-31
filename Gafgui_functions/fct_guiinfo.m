% --------------------------------------------------------------------
function str = fct_guiinfo(handles)

%%%%%%%%%%%%
str = sprintf('%s\n',handles.ifilename);
if ~strcmp(handles.cfilename,'')
    str = sprintf('%s[CALCURVE: %s]\n',str,handles.cfilename);
else
    str = sprintf('%s[CALCURVE: N/D]\n',str);
end
if handles.data==1
    dumb='SIGNAL';
elseif handles.data==2
    dumb='Raw OD';
elseif handles.data==3
    dumb='Net OD';
elseif handles.data==4
    dumb='DOSE';
end
str=sprintf('%s   [Data = %s]',str,dumb);
%%%%%%%%%%%%
str=sprintf('%s   [Filter = %.1f mm]',str,handles.filter*10);
if handles.data==1
    str=sprintf('%s   [Saturation = N/D]',str);
    %elseif handles.data==2
    %    str=sprintf('%s   [Saturation = N/D]',str);
elseif (handles.data==2)||(handles.data==3)
    if handles.SATMAX==0
        str=sprintf('%s   [Saturation = N/D]',str);
    else
        str=sprintf('%s   [Saturation = %.3f - %.3f]',str,handles.SATMIN,handles.SATMAX);
    end
else
    if handles.SATMAX==0
        str=sprintf('%s   [Saturation = N/D]',str);
    else
        str=sprintf('%s   [Saturation = %.1f - %.1f]',str,handles.SATMIN,handles.SATMAX);
    end
end
if max(size(handles.DOSEREF))==0
    str=sprintf('%s   [Reference dose = N/D]\n',str);
else
    str=sprintf('%s   [Reference dose = 1 value]\n',str);
end
if handles.channel==4
    dumb='Multi';
elseif handles.channel==5
    dumb='RGB'
else
    dumb = fct_channeltocolor(handles.channel);
end

str=sprintf('%s   [Channel = %s]',str,dumb);
%%%%%%%%%%%%
str=sprintf('%s   [Bits = %d]',str,handles.BITSini);
%%%%%%%%%%%%
str=sprintf('%s   [Res. = %.3f mm]',str,handles.DELTA*10);
if max(size(handles.BCKGRND))==0
    str=sprintf('%s   [Bckgrnd OD = N/D]',str);
else
    str=sprintf('%s   [Bckgrnd OD = %.0f value(s)]',str,max(size(handles.BCKGRND)));
end