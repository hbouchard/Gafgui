% --------------------------------------------------------------------
function [v,nkernels,index] = fct_vandermatrixRGBfinal(R,G,B,order,zero,negative)

x = R(:); y = G(:); z = B(:);

index = [];
if negative
    k1 = -order;
else
    k1 = 0;
end
k2 = order;
for i = k1:k2
    for j = k1:k2
        for k = k1:k2
            %for l = k1:k2
            index = cat(1,index,[i j k]);
            %end
        end
    end
end
s = sum(index,2);
kk1 = intersect(find(s<=k2),find(s>0));
kk2 = intersect(find(s>=k1),find(s<0));
k = sort([kk1(:)' kk2(:)']);
index = index(k,:);
s = s(k);

v = [];
if zero
    v = cat(2,v,ones(length(x),1));
end
for j = 1:size(index,1)
    v = cat(2,v,[x.^index(j,1).*y.^index(j,2).*z.^index(j,3)]);
end
nkernels = length(k);