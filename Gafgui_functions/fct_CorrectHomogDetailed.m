function [Imcorr,xk] = fct_CorrectHomogDetailed(Im,corrres,rawimsize,corrdir,xkrange,rrange,grange,brange,pr,pg,pb)
%% Make sure x match to xunique
%Here by convention I use x is the direction of correction
%Here I set the min and max values of xunique given by the homog calibration
flag = 1;
%18 March 2020: We need a test here to assure that the image size is never compromised.
[lin,col,depth] = size(Im);
Imcorr = [];  
xk = [];
if depth~=3
    error('Wrong type of image, should be rgb.');
    flag = 0;
elseif (rawimsize(1)~=lin)||(rawimsize(2)~=col)
    tmp = sprintf('Wrong image size, it should be %d x %d. The resolution should be %f cm.',rawimsize(1),rawimsize(2),corrres);
    error(tmp);
    flag = 0;
end
if flag==1
    xk = xkrange(1):xkrange(2);
    %make image restrained to valid positions    
    if corrdir==1
    %     R = R(k,:);
    %     G = G(k,:);
    %     B = B(k,:);
        R = double(Im(xk,:,1));
        G = double(Im(xk,:,2));
        B = double(Im(xk,:,3));
    elseif corrdir==2
    %     R = R(:,k);
    %     G = G(:,k);
    %     B = B(:,k);
        R = double(Im(:,xk,1));
        G = double(Im(:,xk,2));
        B = double(Im(:,xk,3));    
    end
    %new width
    K = length(xk);
    %% Do correction
    %put everything into column vectors
    [lin,col] = size(R);
    Rcorr = zeros(lin,col);
    Gcorr = zeros(lin,col);
    Bcorr = zeros(lin,col);
    is = [0:65535]';
    %HB March 2020: at the moment there isn't anything to limit R,G,B to
    %the ranges used during corrmat creation - rrange, grange, brange
    h = waitbar(0,'Please wait');
    for k = 1:K
        waitbar(k/K,h);
        %here we don't need to worry about min and max values since the
        %function is linear in signal; so no clipping
        %we can also use the positions within the range since the correction
        %matrix is continuous
        ir = min(65535,max(0,pr(1,k).*is.^2+pr(2,k).*is.^1+pr(3,k).*is.^0));
        ig = min(65535,max(0,pg(1,k).*is.^2+pg(2,k).*is.^1+pg(3,k).*is.^0));
        ib = min(65535,max(0,pb(1,k).*is.^2+pb(2,k).*is.^1+pb(3,k).*is.^0));
        [~,irunique] = unique(ir);
        [~,igunique] = unique(ig);
        [~,ibunique] = unique(ib);
        if corrdir==1
            rtmp = max(min(ir),min(max(ir),R(k,:)));
            gtmp = max(min(ig),min(max(ig),G(k,:)));
            btmp = max(min(ib),min(max(ib),B(k,:)));
        elseif corrdir==2
            rtmp = max(min(ir),min(max(ir),R(:,k)));
            gtmp = max(min(ig),min(max(ig),G(:,k)));
            btmp = max(min(ib),min(max(ib),B(:,k)));        
        end
        rcorr = interp1(ir(irunique),is(irunique),rtmp);
        gcorr = interp1(ig(igunique),is(igunique),gtmp);
        bcorr = interp1(ib(ibunique),is(ibunique),btmp);
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
    
    clear Imcorr;
    
    Imcorr(:,:,1) = uint16(Rcorr);
    Imcorr(:,:,2) = uint16(Gcorr);
    Imcorr(:,:,3) = uint16(Bcorr);
end