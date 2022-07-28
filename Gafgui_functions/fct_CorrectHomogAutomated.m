function [Imcorr,xcorr] = fct_CorrectHomogAutomated(x,I,PR,PG,PB,xunique,corrdir)
%% The model
%the 2D model; observed signal = f(position,corrected signal)
%this is the long way to write it
% F = @(x,s) [s.^0.*x.^0 s.^0.*x.^1 s.^0.*x.^2 s.^0.*x.^3 ...
%      s.^0.*x.^4 s.^0.*x.^5  ...
%      s.^1.*x.^0 s.^1.*x.^1 s.^1.*x.^2 s.^1.*x.^3 ...
%      s.^1.*x.^4 s.^1.*x.^5 ];
%this is the compact/general way to write it 
% F = @(x,s,n) cat(2,repmat(s(:),1,n+1).^repmat(zeros(1,n+1),length(s(:)),1).*repmat(x(:),1,n+1).^repmat(0:n,length(x(:)),1),...
%                    repmat(s(:),1,n+1).^repmat(ones(1,n+1) ,length(s(:)),1).*repmat(x(:),1,n+1).^repmat(0:n,length(x(:)),1)    );
% orderR = length(PR)/2-1;
% orderG = length(PG)/2-1;
% orderB = length(PB)/2-1;
%28 July 2020: based on Van Battum 2016 and Schoenfeld 2016, I make it symmetric with x
F = @(x,s,n) cat(2,repmat(s(:),1,n/2+1).^repmat(zeros(1,n/2+1),length(s(:)),1).*repmat(x(:),1,n/2+1).^repmat(0:2:n,length(x(:)),1),...
                          repmat(s(:),1,n/2+1).^repmat(ones(1,n/2+1) ,length(s(:)),1).*repmat(x(:),1,n/2+1).^repmat(0:2:n,length(x(:)),1)    );
orderR = (length(PR)/2-1)*2;
orderG = (length(PG)/2-1)*2;
orderB = (length(PB)/2-1)*2;
%%
R = double(I(:,:,1));
G = double(I(:,:,2));
B = double(I(:,:,3));
%% Make sure x is within the range of xunique
%Here by convention I use x is the corrcorrdirection of correction
%Here I set the min and max values of xunique given by the homog calibration
xmin = min(xunique); 
xmax = max(xunique);
%%%%%%%%%%%%%%%%%%%%
%HACK!
%HB 12 July: this is a hack to avoid width limitation
button = questdlg('Do you want to extend the correction to the whole scanner bed?','Cheat','Yes','No','Yes') ;
if strcmp(button,'Yes')
    xmin = -100; xmax = 100;
end
%%%%%%%%%%%%%%%%%%%%
%
k = intersect(find(x>=xmin),find(x<=xmax));
if corrdir==1
    R = R(k,:);
    G = G(k,:);
    B = B(k,:);
elseif corrdir==2
    R = R(:,k);
    G = G(:,k);
    B = B(:,k);
end
%new width
K = length(k);
xcorr = x(k);
clear k;
%% Do correction
%put everything into column vectors
[lin,col] = size(R);
Rcorr = zeros(lin,col);
Gcorr = zeros(lin,col);
Bcorr = zeros(lin,col);
is = [0:65535]';

h = waitbar(0,'Please wait');
for k = 1:K
    waitbar(k/K,h);
    i = find(xcorr == xcorr(k));
    %here we don't need to worry about min and max values since the
    %function is linear in signal; so no clipping
    %we can also use the positions within the range since the correction
    %matrix is continuous
    ir = F(xcorr(k),is,orderR)*PR;
    ig = F(xcorr(k),is,orderG)*PG;
    ib = F(xcorr(k),is,orderB)*PB;
    if corrdir==1
        rtmp = max(min(ir),min(max(ir),R(i,:)));
        gtmp = max(min(ig),min(max(ig),G(i,:)));
        btmp = max(min(ib),min(max(ib),B(i,:)));
    elseif corrdir==2
        rtmp = max(min(ir),min(max(ir),R(:,i)));
        gtmp = max(min(ig),min(max(ig),G(:,i)));
        btmp = max(min(ib),min(max(ib),B(:,i)));        
    end
    rcorr = min(65535,max(0,interp1(ir,is,rtmp)));
    gcorr = min(65535,max(0,interp1(ig,is,gtmp)));
    bcorr = min(65535,max(0,interp1(ib,is,btmp)));
    if corrdir==1
        Rcorr(k,:) = rcorr(:)';
        Gcorr(k,:) = gcorr(:)';
        Bcorr(k,:) = bcorr(:)';
    elseif corrdir==2
        Rcorr(:,k) = rcorr(:);
        Gcorr(:,k) = gcorr(:);
        Bcorr(:,k) = bcorr(:);     
    end
end
close(h);

Imcorr(:,:,1) = uint16(Rcorr);
Imcorr(:,:,2) = uint16(Gcorr);
Imcorr(:,:,3) = uint16(Bcorr);