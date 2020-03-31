% --------------------------------------------------------------------
function answer = fct_isthereanimage(handles)

answer = 0;
if isfield(handles,'Z')
    if max(size(handles.Z))>0
        answer = 1;
    end
end
