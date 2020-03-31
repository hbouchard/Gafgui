% --------------------------------------------------------------------
function color = fct_MODEtocolor(channel)

if channel==1
    color = [0.8 0.5 0.5];
elseif channel==2
    color = [0.5 0.8 0.5];
elseif channel==3
    color = [0.5 0.5 0.8];
else
    color = [0 0 0];
end
