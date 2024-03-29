function [posx,posy,rmat,bmat,gmat,stackdir] = fct_ReduceBands(posx,posy,rmat,bmat,gmat,nbbands);

minx = min(min(posx{1})); maxx = max(max(posx{1}));
miny = min(min(posy{1})); maxy = max(max(posy{1}));

for i=2:nbbands
    minx = max(min(min(posx{i})),minx); maxx = min(max(max(posx{i})),maxx);
    miny = max(min(min(posy{i})),miny); maxy = min(max(max(posy{i})),maxy);
end
stackdir = questdlg('In which direction do you stack bands?','Direction','Horizontal','Vertical','Horizontal') ;
for i=1:nbbands
    if strcmp(stackdir,'Vertical')
        px = mean(posx{i},1); 
        k = intersect(find(px>=minx),find(px<=maxx));
        A = posx{i}; B = posy{i}; C = rmat{i}; D = gmat{i}; E = bmat{i};
        posx{i} = A(:,k); posy{i} = B(:,k); rmat{i} = C(:,k); gmat{i} = D(:,k); bmat{i} = E(:,k);
    else
        py = mean(posy{i},2); 
        k = intersect(find(py>=miny),find(py<=maxy));
        A = posx{i}; B = posy{i}; C = rmat{i}; D = gmat{i}; E = bmat{i};
        posx{i} = A(k,:); posy{i} = B(k,:); rmat{i} = C(k,:); gmat{i} = D(k,:); bmat{i} = E(k,:);        
    end    
end
clear A B C D E;

R = []; G = []; B = []; pos = [];
for i=1:nbbands
    if strcmp(stackdir,'Vertical')
        R = cat(1,R,rmat{i}); G= cat(1,G,gmat{i}); B = cat(1,B,bmat{i});
        pos = cat(1,pos,posx{i});
    else
        R = cat(2,R,rmat{i}); G= cat(2,G,gmat{i}); B = cat(2,B,bmat{i});  
        pos = cat(2,pos,posy{i});
    end
end
[m,n] = size(R); 
if strcmp(stackdir,'Vertical')
    posx1 = mean(pos,1);
    posy1 = (0:(m-1)-(m-1)/2);
else
    posx1 = (0:(n-1)-(n-1)/2);
    posy1 = mean(pos,2);
end
I(:,:,1) = R; I(:,:,2) = G; I(:,:,3) = B;

figure; imagesc(posx1,posy1,I); set(gca,'dataaspectratio',[1 1 1]);
