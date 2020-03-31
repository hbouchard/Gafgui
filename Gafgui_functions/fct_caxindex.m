

% --------------------------------------------------------------------
function [i,j,k]=fct_caxindex(y)

if size(y,1)>size(y,2)
    y=y';
end
if size(y,1)>1
    error('Function only valid for non-empty arrays.');
end
if size(y,2)<3
    error('Profile too small.');
end
dleft=zeros(1,size(y,2));
dright=zeros(1,size(y,2));

for i=2:(size(y,2)-1)
    halfval=y(i)/2;
    yright(1:size(y,2))=y(i);
    yleft(1:size(y,2))=y(i);
    yright(i:size(y,2))=y(i:size(y,2));
    yleft(1:i)=y(1:i);
    yright=yright-halfval;
    yleft=yleft-halfval;
    
    if min(yright)<0
        iright=max(i,min(find(abs(yright)==min(abs(yright)))));
    else
        iright=i;
    end
    if min(yleft)<0
        ileft =min(i,max(find(abs(yleft)==min(abs(yleft)))));
    else
        ileft=i;
    end
    
    dleft(i)=i-ileft;
    dright(i)=iright-i;
    
    if dleft(i)==0
        dleft(i)=-size(y,2);
    end
    if dright(i)==0
        dright(i)=-2*size(y,2);
    end
    
end
diff=abs(dleft-dright);
dleft=max(dleft,0);
dright=max(dright,0);
diff(1)=2*size(y,2);
diff(size(y,2))=diff(1);
if min(diff)>(size(y,2)-1)
    i=1;
    j=floor((size(y,2)-1)/2);
    k=size(y,2);
else
    j=min(find(diff==min(diff)));
    i=j-dleft(j);
    k=dright(j)+j;
end
