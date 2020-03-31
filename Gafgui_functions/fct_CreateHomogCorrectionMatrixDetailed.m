function handles = fct_CreateHomogCorrectionMatrixDetailed(handles);

% clc;
% clear all;
% close all;
% 
% fct_AddGafguiFctPath();
% handles = fct_initGafgui();
%%

button = questdlg('Do you want to start from the beginning?','Shortcut','Yes','No','Yes') ;
figure(handles.H);

if strcmp(button,'No')
    flag = 1;
elseif strcmp(button,'Yes')
    flag = 2;
else
    flag = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%Get data
if flag == 2%not the shortcut
    
    [posx,posy,rmat,bmat,gmat,res,nbbands,stackdir,bits,posxk,posyk,rawimsize] = fct_GetBandsAuto(handles);
    
    if nbbands==0
        flag = -1;
    elseif nbbands<=2
        errordlg('The number of bands must be at least 3.');
        flag = -1;
    else
        %put band in increasing OD order
        mu = zeros(nbbands,1);
        for i=1:nbbands
           mu(i) = mean(mean(rmat{i})); 
        end
        [mu,index] = sort(mu);
        index = flipud(index(:));
        for i=1:nbbands
           rmat_{i} = rmat{index(i)};
           gmat_{i} = gmat{index(i)};
           bmat_{i} = bmat{index(i)};
           posx_{i} = posx{index(i)};
           posy_{i} = posy{index(i)};
           posxk_{i} = posxk{index(i)};
           posyk_{i} = posyk{index(i)};
        end
        rmat = rmat_; clear rmat_;
        gmat = gmat_; clear gmat_;
        bmat = bmat_; clear bmat_;
        posx = posx_; clear posx_;
        posy = posy_; clear posy_;        
        posxk = posxk_; clear posxk_;
        posyk = posyk_; clear posyk_;
        %%%
        button = questdlg('Do you want to save this step?','Shortcut','Yes','No','Yes') ;
        if strcmp(button,'Yes')
            %ici ajouter un while au cas ou cela ne marche pas
            [ofilename,opathname]=uiputfile({'*.mat'},'Work to save');
            if ofilename==0
                flag == -1;
            else
                filename = fct_makecleanfilename(opathname,ofilename);
                save(filename,'posx','posy','posxk','posyk','rmat','bmat','gmat','bits','res','rawimsize');
            end
        end
    end
elseif flag == 1
    [ifilename,ipathname]=uigetfile({'*.mat'},'Work to import');
    if ~strcmp(class(ifilename),'double')
        filename = fct_makecleanfilename(ipathname,ifilename);
        load (filename);
    else
        flag = -1;
    end
end

%%
if flag~=-1
    %% Reduce size to identical
    lin = zeros(numel(rmat),1);
    col = zeros(numel(rmat),1);
    for i = 1:numel(rmat)
        [lin(i),col(i)] = size(posx{i});
    end
    if std(lin)==0
        stackdir = 'horizontal';
    elseif std(col)==0
        stackdir = 'vertical';
    else
        error('Something went wrong. Dimensions of bands is onconsistent.');
    end
    %this works for both directions
    lin = min(lin);
    col = min(col);
    for i = 1:numel(rmat)
        posx{i} = posx{i}(1:lin,1:col);
        posy{i} = posy{i}(1:lin,1:col);
        posxk{i} = posxk{i}(1:lin,1:col);
        posyk{i} = posyk{i}(1:lin,1:col);
        rmat{i} = rmat{i}(1:lin,1:col);
        gmat{i} = gmat{i}(1:lin,1:col);
        bmat{i} = bmat{i}(1:lin,1:col);
    end
    %% Detailed correction: assemble matrix
    %positions
    xvect = []; 
    %reference signal values
    r0vect = []; g0vect = []; b0vect = [];
    xvect = []; xkvect = []; 
    rvect = []; gvect = []; bvect = [];
    frcorr = []; fgcorr = []; fbcorr = [];
    R = []; G = []; B = []; X = []; Xk = [];
    R0 = []; G0 = []; B0 = [];
    for i = 1:numel(rmat)
        if strcmp(stackdir,'vertical')
            [lin,col] = size(posx{i});
            x = reshape(posx{i},lin*col,1);
            xk = reshape(posxk{i},lin*col,1);
            dir = 1;
            corrdir = 2;
        elseif strcmp(stackdir,'horizontal')
            [lin,col] = size(posy{i});
            x = reshape(posy{i},lin*col,1);
            xk = reshape(posyk{i},lin*col,1);
            dir = 2;
            corrdir = 1;
        end
        r = reshape(double(rmat{i}),lin*col,1);
        g = reshape(double(gmat{i}),lin*col,1);
        b = reshape(double(bmat{i}),lin*col,1);
        %the central 1 cm is the region we consider to be the reference/true signal
        %k = intersect(find(x>=-0.5),find(x<=0.5));
        k = intersect(find(x>=-2.0),find(x<=2.0));
        r0 = 0*r + mean(r(k));
        g0 = 0*g + mean(g(k));
        b0 = 0*b + mean(b(k));
        %
        xkvect = cat(1,xkvect,xk);
        xvect = cat(1,xvect,x);
        rvect = cat(1,rvect,r);     
        gvect = cat(1,gvect,g);       
        bvect = cat(1,bvect,b);
        r0vect = cat(1,r0vect,r0);
        g0vect = cat(1,g0vect,g0);
        b0vect = cat(1,b0vect,b0);
        %
        Xk = cat(dir,Xk,reshape(xk,lin,col));
        X = cat(dir,X,reshape(x,lin,col));
        R = cat(dir,R,reshape(r,lin,col));
        G = cat(dir,G,reshape(g,lin,col));
        B = cat(dir,B,reshape(b,lin,col));
        R0 = cat(dir,R0,reshape(r,lin,col)*0+mean(r0));
        G0 = cat(dir,G0,reshape(g,lin,col)*0+mean(g0));
        B0 = cat(dir,B0,reshape(b,lin,col)*0+mean(b0));
    end
    %% Detailed correction: off-axis each position gets its own polynomial 
    xunique = sort(unique(xvect));
    xkunique = sort(unique(xkvect));
    pr = []; pg = []; pb = [];
    %just a random one to show
    krand = max(1,floor(rand(1,1)*length(xkunique)));
    h = waitbar(0,'Patience');
    for k = 1:length(xkunique)
        waitbar(k/length(xkunique),h);
        if strcmp(stackdir,'vertical')
            r = R(:,k); r0 = R0(:,k);
            g = G(:,k); g0 = G0(:,k);
            b = B(:,k); b0 = B0(:,k);
        elseif strcmp(stackdir,'horizontal')
            r = R(k,:); r0 = R0(k,:);
            g = G(k,:); g0 = G0(k,:);
            b = B(k,:); b0 = B0(k,:);       
        end 
        
        %The fit
        N = 1;%We limit the fit to N = 1 (linear) or 2 (quadratic)
        %
        ptmp0 = [0;1;0];
        ptmp1 = polyfit(r0,r,N)';
        ptmp2 = polyfit(g0,g,N)';
        ptmp3 = polyfit(b0,b,N)';
        if N==1%make sure it is of quadratic form
            ptmp1 = cat(1,0,ptmp1);
            ptmp2 = cat(1,0,ptmp2);
            ptmp3 = cat(1,0,ptmp3);
        end
        if (xunique(k)>=-2.0)&&(xunique(k)<=2.0);
            pr = cat(2,pr,ptmp0);
            pg = cat(2,pg,ptmp0);
            pb = cat(2,pb,ptmp0);
%         elseif (xunique(k)>=-2.0)&&(xunique(k)<=2.0);
%             dist = min([abs(xunique(k)+2.0) abs(xunique(k)-2.0)]);
%             w = (1-dist)/2;
%             tmp0 = [0;1;0];
%             pr = cat(2,pr,w*ptmp0+(1-w)*ptmp1);
%             pg = cat(2,pg,w*ptmp0+(1-w)*ptmp2);
%             pb = cat(2,pb,w*ptmp0+(1-w)*ptmp3);
        else
            pr = cat(2,pr,ptmp1);
            pg = cat(2,pg,ptmp2);
            pb = cat(2,pb,ptmp3);
        end
        %A random show of the fit
        if k==krand
            sr = [min(R(:)):100:max(R(:))];
            sg = [min(G(:)):100:max(G(:))];
            sb = [min(B(:)):100:max(B(:))];
            if 1
                figure;
                hold on;
                plot(r0,r,'.r',sr,min(65535,max(0,polyval(polyfit(r0,r,N),sr))),'-r');
                plot(g0,g,'.g',sg,min(65535,max(0,polyval(polyfit(g0,g,N),sg))),'-g');
                plot(b0,b,'.b',sb,min(65535,max(0,polyval(polyfit(b0,b,N),sb))),'-b');
                hold off;
            else
                Fr = [r0(:).^3 r0(:).^2 r0(:).^1]; pr = Fr\r(:);
                Fg = [g0(:).^3 g0(:).^2 g0(:).^1]; pg = Fg\g(:);
                Fb = [b0(:).^3 b0(:).^2 b0(:).^1]; pb = Fb\b(:);
                Fr = [sr(:).^3 sr(:).^2 sr(:).^1];
                Fg = [sg(:).^3 sg(:).^2 sg(:).^1];
                Fb = [sb(:).^3 sb(:).^2 sb(:).^1];
                figure;
                hold on;
                plot(r0,r,'.r',sr,min(65535,max(0,Fr*pr)),'-r');
                plot(g0,g,'.g',sg,min(65535,max(0,Fg*pg)),'-g');
                plot(b0,b,'.b',sb,min(65535,max(0,Fb*pb)),'-b');
                hold off;
            end
        end
    end
    close(h);
 %% Vizualise correction map   
    sr = [min(R(:)):100:max(R(:))]';
    sg = [min(G(:)):100:max(G(:))]';
    sb = [min(B(:)):100:max(B(:))]';
    Sr = repmat(sr,1,length(xkunique));
    Sg = repmat(sg,1,length(xkunique));
    Sb = repmat(sb,1,length(xkunique));
    Xr = repmat(xunique(:)',length(sr),1);
    Xg = repmat(xunique(:)',length(sg),1);
    Xb = repmat(xunique(:)',length(sb),1);
    PR1 = repmat(pr(1,:),length(sr),1);
    PR2 = repmat(pr(2,:),length(sr),1);
    PR3 = repmat(pr(3,:),length(sr),1);
    PG1 = repmat(pg(1,:),length(sg),1);
    PG2 = repmat(pg(2,:),length(sg),1);
    PG3 = repmat(pg(3,:),length(sg),1); 
    PB1 = repmat(pb(1,:),length(sb),1);
    PB2 = repmat(pb(2,:),length(sb),1);
    PB3 = repmat(pb(3,:),length(sb),1);  
    Rhat = min(65535,max(0,PR1.*Sr.^2+PR2.*Sr.^1+PR3.*Sr.^0));
    Ghat = min(65535,max(0,PG1.*Sg.^2+PG2.*Sg.^1+PG3.*Sg.^0));
    Bhat = min(65535,max(0,PB1.*Sb.^2+PB2.*Sb.^1+PB3.*Sb.^0));

    figure;
    subplot(1,3,1);
    mesh(Xr,Sr,Rhat./Sr);
    title('Red');
    subplot(1,3,2);
    mesh(Xg,Sg,Ghat./Sg);
    title('Green');
    subplot(1,3,3);
    mesh(Xb,Sb,Bhat./Sb);
    title('Blue');
    %%
    Xk = []; R = []; G = []; B = []; 
    for i =1:numel(rmat)
        [lin,col] = size(rmat{i});
        if lin<col %Horizontal correction
            Xk = cat(1,double(posxk{i}),Xk);
            R = cat(1,double(rmat{i}),R);
            G = cat(1,double(gmat{i}),G);
            B = cat(1,double(bmat{i}),B);
        elseif lin>col %Vertical correction
            Xk = cat(2,double(posyk{i}),Xk);
            R = cat(2,double(rmat{i}),R);
            G = cat(2,double(gmat{i}),G);
            B = cat(2,double(bmat{i}),B);
        end
    end
    %
    if corrdir==1
        xk = Xk(:,1);
    elseif corrdir==2
        xk = Xk(1,:);
    end
    Im(:,:,1) = uint16(R);
    Im(:,:,2) = uint16(G);
    Im(:,:,3) = uint16(B);

    %%
    xrange = [min(xunique) max(xunique)];
    xkrange = [min(xkunique) max(xkunique)];
    rrange = [min(R(:)) max(R(:))];
    grange = [min(G(:)) max(G(:))];
    brange = [min(B(:)) max(B(:))];
    %Here this is a special use of the function since Im is already reduced
    [n,m] = size(Im(:,:,1));
    Imcorr = fct_CorrectHomogDetailed(Im,res,[n m],corrdir,[1 (1+max(xkrange)-min(xkrange))],rrange,grange,brange,pr,pg,pb);
    %%
    figure;
    subplot(3,2,1); mesh(Im(:,:,1));     title('Red uncorrected');
    subplot(3,2,2); mesh(Imcorr(:,:,1)); title('Red corrected');
    subplot(3,2,3); mesh(Im(:,:,2));     title('Green uncorrected');
    subplot(3,2,4); mesh(Imcorr(:,:,2)); title('Green corrected');
    subplot(3,2,5); mesh(Im(:,:,3));     title('Blue uncorrected');
    subplot(3,2,6); mesh(Imcorr(:,:,3)); title('Blue corrected');
    %% Check with Gafgui if corrected is fine 
%     tmp = 'SomeName.tif';
%     
%     defaultname = sprintf('%s_homog.tif',tmp(1:(length(tmp)-4)));
%     [ofilename,opathname] = uiputfile({'*.tif'},'Save corrected image',defaultname);
% 
%     if (~strcmp(class(ofilename),'double'))
%         filename = fct_makecleanfilename(opathname,ofilename);
%         imwrite(Imcorr,filename,'tif','Resolution',round(2.54/res));
%     end       
    %% Write correction matrix
    [ofilename,opathname]=uiputfile({'*.hm3'},['Save homogeneity correction']);
    if ~strcmp(class(ofilename),'double')
        fname = fct_makecleanfilename(opathname,ofilename);
        err = fct_WriteCorrmatrixDetailed(fname,res,rawimsize,corrdir,xkrange,rrange,grange,brange,pr,pg,pb);
        if err
            errordlg('File not written');
        else
            [corrres_,rawimsize_,corrdir_,xkrange_,rrange_,grange_,brange_,pr_,pg_,pb_]= fct_ReadCorrmatrixDetailed(fname);
            diff = [   corrres_-res;
                rawimsize_(1)-rawimsize(1);
                rawimsize_(2)-rawimsize(2);
                xkrange_(1) - xkrange(1);
                xkrange_(2) - xkrange(2);
                rrange_(1) - rrange(1);
                rrange_(2) - rrange(2);
                grange_(1) - grange(1);
                grange_(2) - grange(2);
                brange_(1) - brange(1);
                brange_(2) - brange(2);
                max(max(abs(pr_-pr)));
                max(max(abs(pg_-pg)));
                max(max(abs(pb_-pb)));];
            diff
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%