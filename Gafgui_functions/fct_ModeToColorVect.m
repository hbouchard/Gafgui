% --------------------------------------------------------------------
function color = fct_ModeToColorVect(channel)

%https://www.rapidtables.com/web/color/RGB_Color.html

if channel==1 %Red
%     color = [0.8 0.5 0.5]; 
    color = [175 31 31]/255;
elseif channel==2 %Green
%     color = [0.5 0.8 0.5];
    color = [31 175 31]/255;
elseif channel==3 %Blue
%     color = [0.5 0.5 0.8];
    color = [31 31 175]/255;
elseif channel==4 %Multi
%     color = [0.8 0.8 0.5];    
    color = [175 175 31]/255;
else %Unknown
    color = [0.5 0.5 0.5];
end
