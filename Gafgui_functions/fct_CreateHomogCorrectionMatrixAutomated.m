function handles = fct_CreateHomogCorrectionMatrixAutomated(handles);
% 
% clc;
% clear all;
% close all;
% 
% fct_AddGafguiFctPath();
% handles = fct_initGafgui();
%% Comments HB 11 March 2020
% Some important points: 
%
% 1. the choice of the order of the fit can be defined by the user in terms of position 
% dependence only. I recommend using linear only for signal, based both on observation 
% but also not to restrain/clip the image to min/max values used in the
% calibration process
%
% 2. the possibility to correct, no matter what the resolution is, within the  
% position boundaries. At the moment, the correction method force the image
% to have a given value of position.
%
%%
button = questdlg('Do you want to start from the beginning?','Shortcut','Yes','No','Yes') ;
% figure(handles.H);

if strcmp(button,'No')
    flag = 1;
else
    flag = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%Get data
if flag == 2%shortcut
    
    [posx,posy,rmat,bmat,gmat,res,nbbands,stackdir,bits] = fct_GetBandsAuto(handles);
 
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
        end
        rmat = rmat_; clear rmat_;
        gmat = gmat_; clear gmat_;
        bmat = bmat_; clear bmat_;
        posx = posx_; clear posx_;
        posy = posy_; clear posy_;
        %%%
        button = questdlg('Do you want to save this step?','Shortcut','Yes','No','Yes') ;
        if strcmp(button,'Yes')
            %ici ajouter un while au cas ou cela ne marche pas
            [ofilename,opathname]=uiputfile({'*.mat'},'Work to save');
            if ofilename==0
                flag == -1;
            else
                filename = fct_makecleanfilename(opathname,ofilename);
                save(filename,'posx','posy','rmat','bmat','gmat','bits','res');
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
        rmat{i} = rmat{i}(1:lin,1:col);
        gmat{i} = gmat{i}(1:lin,1:col);
        bmat{i} = bmat{i}(1:lin,1:col);
    end
    %% Detailed correction
    %positions
    xvect = []; 
    %reference signal values
    r0vect = []; g0vect = []; b0vect = [];
    rvect = []; gvect = []; bvect = [];
    frcorr = []; fgcorr = []; fbcorr = [];
    R = []; G = []; B = []; X = [];
    R0 = []; G0 = []; B0 = [];
    for i = 1:numel(rmat)
        if strcmp(stackdir,'vertical')
            [lin,col] = size(posx{i});
            x = reshape(posx{i},lin*col,1);
            dir = 1;
        else
            [lin,col] = size(posy{i});
            x = reshape(posy{i},lin*col,1);
            dir = 2;
        end
        r = reshape(double(rmat{i}),lin*col,1);
        g = reshape(double(gmat{i}),lin*col,1);
        b = reshape(double(bmat{i}),lin*col,1);
        %the central 1 cm is the region we consider to be the reference/true signal
        k = intersect(find(x>=-0.5),find(x<=0.5));
        r0 = 0*r + mean(r(k));
        g0 = 0*g + mean(g(k));
        b0 = 0*b + mean(b(k));
        %
        xvect = cat(1,xvect,x);
        rvect = cat(1,rvect,r);     
        gvect = cat(1,gvect,g);       
        bvect = cat(1,bvect,b);
        r0vect = cat(1,r0vect,r0);
        g0vect = cat(1,g0vect,g0);
        b0vect = cat(1,b0vect,b0);
        %
        X = cat(dir,X,reshape(x,lin,col));
        R = cat(dir,R,reshape(r,lin,col));
        G = cat(dir,G,reshape(g,lin,col));
        B = cat(dir,B,reshape(b,lin,col));
        R0 = cat(dir,R0,reshape(r,lin,col)*0+mean(r0));
        G0 = cat(dir,G0,reshape(g,lin,col)*0+mean(g0));
        B0 = cat(dir,B0,reshape(b,lin,col)*0+mean(b0));
    end
    %% Display extracted signal relative to center as a function of positions
    %the 2D model; observed signal = f(position,corrected signal)
    %this is the long way to write it
%     F = @(x,s) [s.^0.*x.^0 s.^0.*x.^1 s.^0.*x.^2 s.^0.*x.^3 ...
%          s.^0.*x.^4 s.^0.*x.^5  ...
%          s.^1.*x.^0 s.^1.*x.^1 s.^1.*x.^2 s.^1.*x.^3 ...
%          s.^1.*x.^4 s.^1.*x.^5 ];
    %this is the compact/general way to write it 
    %F = @(x,s,n) cat(2,repmat(s(:),1,n+1).^repmat(zeros(1,n+1),length(s(:)),1).*repmat(x(:),1,n+1).^repmat(0:n,length(x(:)),1),...
    %                   repmat(s(:),1,n+1).^repmat(ones(1,n+1) ,length(s(:)),1).*repmat(x(:),1,n+1).^repmat(0:n,length(x(:)),1)    );
    %28 July 2020: based on Van Battum 2016 and Schoenfeld 2016, I make it symmetric with x
    F = @(x,s,n) cat(2,repmat(s(:),1,n/2+1).^repmat(zeros(1,n/2+1),length(s(:)),1).*repmat(x(:),1,n/2+1).^repmat(0:2:n,length(x(:)),1),...
                          repmat(s(:),1,n/2+1).^repmat(ones(1,n/2+1) ,length(s(:)),1).*repmat(x(:),1,n/2+1).^repmat(0:2:n,length(x(:)),1)    );
                   
    %apply LS fit
    [lin,col] = size(R); 
    x = reshape(X,lin*col,1);
    r0 = reshape(R0,lin*col,1);
    g0 = reshape(G0,lin*col,1);
    b0 = reshape(B0,lin*col,1);
    r = reshape(R,lin*col,1);
    g = reshape(G,lin*col,1);
    b = reshape(B,lin*col,1);
    %fitted coefficients
    unconfirmed = 1;
    while unconfirmed
        % HB 11 March 2020: here the user can choose    
        orderR = -1; orderG = -1; orderB = -1;
        flag = 1;
        clear tmp;
%         while (orderR<1||orderR>6)&&(orderG<1||orderG>6)&&(orderB<1||orderB>6)&&flag==1
%             tmp{1} = 'Order of fit in position for RED (1-5)';
%             tmp{2} = 'Order of fit in position for GREEN (1-5)';
%             tmp{3} = 'Order of fit in position for BLUE (1-5)';
%             answer = inputdlg(tmp,'Input',[1 1 1],{'5','4','4'});
        %28 July 2020: I make sure it is symmetric
        while (orderR<0||orderR>6)||(orderG<0||orderG>6)||(orderB<0||orderB>6) || (mod(orderR,2)~=0||mod(orderG,2)~=0||mod(orderB,2)~=0) && (flag==1)
            tmp{1} = 'Order of fit in position for RED (0,2,4 or 6)';
            tmp{2} = 'Order of fit in position for GREEN (0,2,4 or 6)';
            tmp{3} = 'Order of fit in position for BLUE (0,2,4 or 6)';
            answer = inputdlg(tmp,'Input',[1 1 1],{'4','2','2'});
            if numel(answer)==0
                flag = 0;
                unconfirmed = 0;
            else
                orderR = str2num(answer{1});
                orderG = str2num(answer{2});
                orderB = str2num(answer{3});
            end
        end
        clear tmp;
        %%
        if flag
            PR = F(x,r0,orderR)\r;
            PG = F(x,g0,orderG)\g;
            PB = F(x,b0,orderB)\b;
%           Commented on 28 July 2020
%             %what if I put all colors into 1 fit
%             if 0%does seem like a good idea...
%                 P = F(cat(1,x,cat(1,x,x)),cat(1,r0,cat(1,g0,b0)),orderR)\cat(1,r,cat(1,g,b));
%                 PR = P;
%                 PG = P;
%                 PB = P;
%             %what if I put green and blue into 1 fit
%             elseif 0
%                 P = F(cat(1,x,x),cat(1,g0,b0),orderG)\cat(1,g,b);
%                 PG = P;
%                 PB = P;            
%             end
            %interpolation
            xunique = sort(unique(xvect));
            %here this needs to be fixed
            rlims = [min(R(:)) max(R(:))];
            glims = [min(G(:)) max(G(:))];
            blims = [min(B(:)) max(B(:))];
            r0lims = [min(R0(:)) max(R0(:))];
            g0lims = [min(G0(:)) max(G0(:))];
            b0lims = [min(B0(:)) max(B0(:))];    
            %%%
            [IXr,ISr] = meshgrid(xunique,rlims(1):100:rlims(2));
            [IXg,ISg] = meshgrid(xunique,glims(1):100:glims(2));
            [IXb,ISb] = meshgrid(xunique,blims(1):100:blims(2));
            [linr,colr] = size(IXr); 
            [ling,colg] = size(IXg); 
            [linb,colb] = size(IXb); 
            ixr = reshape(IXr,linr*colr,1);
            ixg = reshape(IXg,ling*colg,1);
            ixb = reshape(IXb,linb*colb,1);
            isr = reshape(ISr,linr*colr,1);
            isg = reshape(ISg,ling*colg,1);
            isb = reshape(ISb,linb*colb,1);
            %predicted correction factors for interpolation
            Fr = reshape(F(ixr,isr,orderR)*PR,linr,colr)./ISr;
            Fg = reshape(F(ixg,isg,orderG)*PG,ling,colg)./ISg;
            Fb = reshape(F(ixb,isb,orderB)*PB,linb,colb)./ISb;
            %Display of raw signal inference
            figure;
            subplot(2,3,1);
            hold on;
            plot3(reshape(X,lin*col,1),reshape(R0,lin*col,1),reshape(R,lin*col,1)./reshape(R0,lin*col,1),'.r','markersize',1);
            mesh(IXr,ISr,Fr); colormap('gray')
            hold off;
            view(-15,15);
            title('Red');
            subplot(2,3,2);
            hold on;
            plot3(reshape(X,lin*col,1),reshape(G0,lin*col,1),reshape(G,lin*col,1)./reshape(G0,lin*col,1),'.g','markersize',1);
            mesh(IXg,ISg,Fg); colormap('gray')
            hold off;
            view(-15,15); 
            title('Green');
            subplot(2,3,3);
            hold on;
            plot3(reshape(X,lin*col,1),reshape(B0,lin*col,1),reshape(B,lin*col,1)./reshape(B0,lin*col,1),'.b','markersize',1);
            mesh(IXb,ISb,Fb); colormap('gray')
            hold off;
            view(-15,15);
            title('Blue');
            %predicted correction factors for interpolation
            [IXr,ISr] = meshgrid(xunique,rlims(1):100:rlims(2));
            [IXg,ISg] = meshgrid(xunique,glims(1):100:glims(2));
            [IXb,ISb] = meshgrid(xunique,blims(1):100:blims(2));
            [linr,colr] = size(IXr); 
            [ling,colg] = size(IXg); 
            [linb,colb] = size(IXb); 
            ixr = reshape(IXr,linr*colr,1);
            ixg = reshape(IXg,ling*colg,1);
            ixb = reshape(IXb,linb*colb,1);
            isr = reshape(ISr,linr*colr,1);
            isg = reshape(ISg,ling*colg,1);
            isb = reshape(ISb,linb*colb,1);
            %predicted correction factors for interpolation
            Fr = reshape(F(ixr,isr,orderR)*PR,linr,colr)./ISr;
            Fg = reshape(F(ixg,isg,orderG)*PG,ling,colg)./ISg;
            Fb = reshape(F(ixb,isb,orderB)*PB,linb,colb)./ISb;
            %Display of raw signal inference
            subplot(2,3,4);
            mesh(IXr,ISr,Fr); colormap('gray')
            view(-15,15);
            title('Red');
            subplot(2,3,5);
            mesh(IXg,ISg,Fg); colormap('gray')
            view(-15,15); 
            title('Green');
            subplot(2,3,6);
            mesh(IXb,ISb,Fb); colormap('gray')
            view(-15,15);
            title('Blue');
%             %28 July 2020: this is a test to make sure it is consisten with
%             %our new paper
%             Ar = 0; Br = 0; 
%             for i=1:length(PR)/2
%                 Ar = Ar + PR(i)*IXr.^(2*(i-1));
%                 Br = Br + PR(length(PR)/2+i)*IXr.^(2*(i-1));
%             end
%             Ag = 0; Bg = 0;
%             for i=1:length(PG)/2
%                 Ag = Ag + PG(i)*IXg.^(2*(i-1));
%                 Bg = Bg + PG(length(PG)/2+i)*IXg.^(2*(i-1));
%             end            
%             Ab = 0; Bb = 0;
%             for i=1:length(PB)/2
%                 Ab = Ab + PB(i)*IXb.^(2*(i-1));
%                 Bb = Bb + PB(length(PB)/2+i)*IXb.^(2*(i-1));
%             end            
%             Ftestr = ISr./((ISr-Ar)./Br);
%             Ftestg = ISg./((ISg-Ag)./Bg);
%             Ftestb = ISb./((ISb-Ab)./Bb);
%             figure;
%             subplot(2,3,1);mesh(IXr,ISr,Ftestr);colormap('jet');view(-15,15);title('Red');
%             subplot(2,3,2);mesh(IXg,ISg,Ftestg);colormap('jet');view(-15,15);title('Green');
%             subplot(2,3,3);mesh(IXb,ISb,Ftestb);colormap('jet');view(-15,15);title('Blue');  
%             subplot(2,3,4);mesh(IXr,ISr,Fr);colormap('jet');view(-15,15);title('Red');
%             subplot(2,3,5);mesh(IXg,ISg,Fg);colormap('jet');view(-15,15);title('Green');
%             subplot(2,3,6);mesh(IXb,ISb,Fb);colormap('jet');view(-15,15);title('Blue');            
            %%
            answer = questdlg('Would you like to keep this fit?', ...
                              'Fit', 'Yes','No','Cancel','Yes');
            % Handle response
            switch answer
                case 'Yes'
                    unconfirmed = 0;
                case 'No'
                    unconfirmed = 1;
                    close;
                case 'Cancel'
                    unconfirmed = 0;    
                    flag = 0;
                    close;
            end
        end
    end
%     %% HB July 2022: testing if it is worth correcting OD inseatd of signal based on van Battum a Schoenfeld
%     itest = 3;
%     ODR = double(rmat{itest});
%     ODR = log10(65535./ODR(:));
%     i = intersect(find(x<1),find(x>-1));
%     xi = mean(ODR(i))./ODR;
%     if dir==2
%         x = posy{itest}; 
%     elseif dir==1
%         x = posx{itest}; 
%     end
%     x = x(:);
%     figure;
%     plot(x,xi,'.'); 
%     xlabel('Position (cm)');
%     ylabel('Signal');
%     %x = l*sqrt(n^2*(1-xi.^2)./(1-n^2.*(1-xi.^2)))+0.5*t*sqrt(1-xi.^2)./xi;
%     hold;
%     xhat = @(xi,l,n,t) l*sqrt(n^2*(1-xi.^2)./(1-n^2.*(1-xi.^2)))+0.5*t*sqrt(1-xi.^2)./xi;
%     ixi = [0.9:0.001:1]';
%     n = 1.2;
%     l = 35.6;
%     t = 50e-4;
%     ixhat = -xhat(ixi,l,n,t);
%     ixhat = cat(1,ixhat,xhat(flipud(ixi(1:end-1)),l,n,t));
%     ixi = cat(1,ixi,flipud(ixi(1:end-1)));
%     plot(ixhat,ixi,'.r');
%     xlim([min(x(:)) max(x(:))]);
% 
%     figure;
%     plot(x(:),interp1(ixhat,ixi,x(:)).*ODR(:),'.k')
%     
%     zzz
    %%
    if flag
        clear r g b x r0 g0 b0;
        %% do correction from R,G,B and X
        %Here we must be careful; stackdir is the orthogonal direction
        %of the correction direction, i.e., if stackdir is horizontal, then
        %correction is vertical and vice-versa. So we do this:
        corrdir = 3-dir; %this means 2->1 and 1->2
        %Here by convention I use x is the direction of correction
        if corrdir==1
            x = X(:,1);
        elseif corrdir==2
            x = X(1,:);   
        end
        Im(:,:,1) = uint16(R);
        Im(:,:,2) = uint16(G);
        Im(:,:,3) = uint16(B);
        
        [Imcorr,xcorr] = fct_CorrectHomogAutomated(x,Im,PR,PG,PB,xunique,corrdir);    
        %% Display
%       28 July 2020: I'm not sure what I am getting from this        
%         figure;
%         subplot(1,2,1);
%         imagesc(R); colormap('gray'); %impixelinfo;
%         subplot(1,2,2);
%         imagesc(Imcorr(:,:,1)); colormap('gray'); %impixelinfo;
% 
%         figure;
%         subplot(1,2,1);
%         imagesc(G); colormap('gray'); %impixelinfo;
%         subplot(1,2,2);
%         imagesc(Imcorr(:,:,2)); colormap('gray'); %impixelinfo;
% 
%         figure;
%         subplot(1,2,1);
%         imagesc(B); colormap('gray'); %impixelinfo;
%         subplot(1,2,2);
%         imagesc(Imcorr(:,:,3)); colormap('gray'); %impixelinfo;
        %% Write correction matrix
        [ofilename,opathname]=uiputfile({'*.hm3'},['Save homogeneity correction']);
        if ~strcmp(class(ofilename),'double')
            fname = fct_makecleanfilename(opathname,ofilename);
            err = fct_WriteHomogCorrAutomated(fname,PR,PG,PB,xunique,rlims,glims,blims);
            if err
                errordlg('File not written');
            elseif 0
                [PR_check,PG_check,PB_check,x_check,rlims,glims,blims]  = fct_ReadHomogCorrAutomated(fname);
                sum(abs(PR-PR_check))
                sum(abs(PR-PR_check))
                sum(abs(PR-PR_check))
                sum(abs(x-x_check))
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%