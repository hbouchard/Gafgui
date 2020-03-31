function  [] = fct_AddGafguiFctPath();

% [status, result] = system('ls');
%
% if status==0
%    c = '/'; % unix
% else
%    c = '\'; % dos
% end

cpu = computer;
if strcmp(cpu,'PCWIN')||strcmp(cpu,'PCWIN64')
    c = '\'; % dos
else
    c = '/'; % mac/linux/unix
    %In the case someone uses a SUN, I don't know if this is correct
end

path = cd;

if (path(max(size(path,1),size(path,2)))==c)
    path = [path 'Gafgui_functions'];
else
    path = [path c 'Gafgui_functions'];
end
addpath(path);
