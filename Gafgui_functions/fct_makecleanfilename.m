% --------------------------------------------------------------------
function name = fct_makecleanfilename(path,filename)

% [status, result] = system('ls');
%
% if status==0
%    c = '/'; % unix
% else
%    c = '\'; % dos
% end
if ischar(filename)&&ischar(path)
    cpu = computer;
    if strcmp(cpu,'PCWIN')||strcmp(cpu,'PCWIN64')
        c = '\'; % dos
    else
        c = '/'; % mac/linux/unix
        %In the case someone uses a SUN, I don't know if this is correct
    end
    
    if (path(max(size(path,1),size(path,2)))==c)
        name = [path filename];
    else
        name = [path c filename];
    end
else
    name = [];
end