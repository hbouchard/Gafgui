% --------------------------------------------------------------------
function color = fct_channeltocolor(channel)

% Return a string corresponding to the code used for the channel
% or the ratio of channels

if channel == 1
    color = 'Red';
elseif channel == 2
    color = 'Green';
elseif channel == 3
    color = 'Blue';
elseif channel == 4
    color = 'Multi';
end