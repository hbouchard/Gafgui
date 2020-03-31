% --------------------------------------------------------------------
function channel = fct_colortochannel(color)

% Returns a numerical value for the channel used

if strcmp(color,'Red')
    channel = 1;
elseif strcmp(color,'Green')
    channel = 2;
elseif strcmp(color,'Blue')
    channel = 3;
elseif strcmp(color,'Multi')
    channel = 4;
else channel = 0;
end