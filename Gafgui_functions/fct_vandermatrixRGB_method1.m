% --------------------------------------------------------------------
function [v,nkernels,index] = fct_vandermatrixRGB_method1(R,G,B,order)

x = R(:); y = G(:); z = B(:);

index = [];
for i = 0:order
    for j = 0:order
        for k = 0:order
            index = cat(1,index,[i j k]);
        end
    end
end

index = index(intersect(find(sum(index,2)<=order),find(sum(index,2)>0)),:);
index = cat(1,[0 0 0],index);
v = [];
for j = 1:size(index,1)
    v = cat(2,v,[x.^index(j,1).*y.^index(j,2).*z.^index(j,3)]);
end

nkernels = size(v,2);