% --------------------------------------------------------------------
function h = fct_display(z,del)

x = [-(size(z,2)-1)/2 (size(z,2)-1)/2].*del;
y = [-(size(z,1)-1)/2 (size(z,1)-1)/2].*del;
imagesc(x,y,z);
if length(size(z))==2
    %we check if the data is in the range of OD values
    val = log10(65535/0.5);
    tmp = z(:);
    tmp = tmp(intersect(find(~isnan(tmp)),find(~isinf(tmp))));
    if max(tmp)<val%ODs
        c = colormap('gray');
        colormap(flipud(c));
    else%signal
        colormap('gray');
    end
end
set(gca,'DataAspectRatio',[1 1 1]);
set(gcf,'Units','pixels');
impixelinfo;
h = gcf;
