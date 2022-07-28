function handles = fct_CreateMultichannelCorrection(handles)
% % 
% clc;
% clear all;
% close all;
% %
% addpath('/Users/Hugo/Documents/MATLAB/work/Gafgui_functions');
% %
% fct_AddGafguiFctPath();
% handles = fct_initGafgui();
% close;

button = questdlg('Do you want to start from the beginning?','Shortcut','Yes','No','Yes') ;
% figure(handles.H);

if strcmp(button,'No')
    flag = 1;%shortcut
else
    flag = 2;%data acquisition
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%Shortcut
if flag == 1
    
    [ifilename,ipathname]=uigetfile({'*.mat'},'Work to import');
    if ~strcmp(class(ifilename),'double')
        filename = fct_makecleanfilename(ipathname,ifilename);
        load (filename);
    else
        flag = -1;
    end
end
%Data acquisition
if flag == 2%data acquisition of scanned films
    
    [posx,posy,rmat,bmat,gmat,res,nbbands,stackdir,bits] = fct_GetBandsAuto(handles);
    %[posx,posy,rmat,bmat,gmat,res,nbbands,stackdir,bits,posxk,posyk,rawimsize] = fct_GetBandsAuto(handles);

    if nbbands==0
        flag = -1;
    elseif nbbands<=2
        errordlg('The number of bands must be at least 3.');
        flag = -1;       
    end
end
%Save environment to continue later from here
if flag==2
    %%%
    button = questdlg('Do you want to save this step?','Shortcut','Yes','No','Yes') ;
    if strcmp(button,'Yes')
        [ofilename,opathname]=uiputfile({'*.mat'},'Work to save');
        if ofilename==0
            flag == -1;
        else
            filename = fct_makecleanfilename(opathname,ofilename);
            save(filename,'posx','posy','rmat','bmat','gmat','stackdir','nbbands');
        end
    end
end
%% We go forward if the user selected to film or mat file to work with
if flag>0
    %% Do PCA
    R = []; G = []; B = []; 
    for i= 1 :nbbands
       if strcmp(stackdir,'Vertical')
            k = 1;
       elseif strcmp(stackdir,'Horizontal')
            k = 2;   
       end
       R = cat(k,R,rmat{i});
       G = cat(k,G,gmat{i});
       B = cat(k,B,bmat{i});                           
    end
    [lin,col] = size(R);
    %signal correction range
    rrange = [min(R(:)) max(R(:))];
    grange = [min(G(:)) max(G(:))];
    brange = [min(B(:)) max(B(:))];
    %convert to OD
    R = reshape(-log10(double(R)/65535),lin*col,1);
    G = reshape(-log10(double(G)/65535),lin*col,1);
    B = reshape(-log10(double(B)/65535),lin*col,1);
    [~,P,~] = fct_GetPCA([R G B],0);
    Rot = P;
    mupc = [mean(R) mean(G) mean(B)]*Rot';

    %The implementation of PCA works this way: % [X Y Z] = [R G B]*Rot';
    %eigenvectors
    E = [1 0 0;
         0 1 0;
         0 0 1;]*Rot';
    v1 = E(:,1); 
    v2 = E(:,2); 
    v3 = E(:,3); 
    %At this point I am pretty sure the first vector is pointing out from origin 
    %into the first octant, unless the data is very badly dependent on dose - 
    %think about it: all OD are positive and increase with increasing dose, so
    %the firs component must be in that direction. I will make sure it is
    %pointing in the right direction. For the other components since they don't 
    %vary much, what matter is their sign
    v1 = sign(mupc(1))*v1; 
    v2 = sign(mupc(2))*v2; 
    v3 = sign(mupc(3))*v3; 
    E = [v1 v2 v3];
    Rot = E';
    tmp = [R(:) G(:) B(:)]*Rot';
    X = tmp(:,1);
    Y = tmp(:,2);
    Z = tmp(:,3);
    clear tmp;
    %
    x0 = mean(R); y0 = mean(G);  z0 = mean(B);
    %
    figure;
    hold on;
    plot3([0 1]*max(R(:)),[0 0]*max(R(:)),[0 0]*max(R(:)),'r','linewidth',2);
    plot3([0 0]*max(G(:)),[0 1]*max(G(:)),[0 0]*max(G(:)),'g','linewidth',2);
    plot3([0 0]*max(B(:)),[0 0]*max(B(:)),[0 1]*max(B(:)),'b','linewidth',2);
    plot3(x0+[-E(1,1) E(1,1)]*max(X(:))/2,y0+[-E(2,1) E(2,1)]*max(X(:))/2,z0+[-E(3,1) E(3,1)]*max(X(:))/2,'color',[1 1 0],'linewidth',2);
    plot3(x0+[-E(1,2) E(1,2)]*max(Y(:))/2,y0+[-E(2,2) E(2,2)]*max(Y(:))/2,z0+[-E(3,2) E(3,2)]*max(Y(:))/2,'color',[0 1 1],'linewidth',2);
    plot3(x0+[-E(1,3) E(1,3)]*max(Z(:))/2,y0+[-E(2,3) E(2,3)]*max(Z(:))/2,z0+[-E(3,3) E(3,3)]*max(Z(:))/2,'color',[1 0 1],'linewidth',2);
    plot3(R(:),G(:),B(:),'.k');
    % x0 = mean(X)*E(1,1); y0 = mean(X)*E(2,1);  z0 = mean(X)*E(3,1);
    % plot3(x0+[0 E(1,2)]*max(Y(:)),...
    %       y0+[0 E(2,2)]*max(Y(:)),...
    %       z0+[0 E(3,2)]*max(Y(:)),'color',[0 1 1],'linewidth',2);
    % plot3(x0+[0 E(1,3)]*max(Z(:)),...
    %       y0+[0 E(2,3)]*max(Z(:)),...
    %       z0+[0 E(3,3)]*max(Z(:)),'color',[1 0 1],'linewidth',2);
    hold off;
    legend('Red','Green','Blue','E#1','E#2','E#3');
    xlabel('x');
    ylabel('y');
    zlabel('z');
    clear R G B X Y Z;

    % imrot = [R G B]*Rot';
    % X = reshape(imrot(:,1), lin,col);
    % Y = reshape(imrot(:,2), lin,col);
    % Z = reshape(imrot(:,3), lin,col);
    % clear imrot;

    r = zeros(nbbands,1);
    g = zeros(nbbands,1);
    b = zeros(nbbands,1);
    x = zeros(nbbands,1);
    y = zeros(nbbands,1);
    z = zeros(nbbands,1);
    for i= 1 :nbbands
        xmat{i} =  -(Rot(1,1)*log10(double(rmat{i})/65535)+Rot(1,2)*log10(double(gmat{i})/65535)+Rot(1,3)*log10(double(bmat{i})/65535));
        ymat{i} =  -(Rot(2,1)*log10(double(rmat{i})/65535)+Rot(2,2)*log10(double(gmat{i})/65535)+Rot(2,3)*log10(double(bmat{i})/65535));
        zmat{i} =  -(Rot(3,1)*log10(double(rmat{i})/65535)+Rot(3,2)*log10(double(gmat{i})/65535)+Rot(3,3)*log10(double(bmat{i})/65535));
        r(i) = mean(mean(-log10(double(rmat{i})/65535)));
        g(i) = mean(mean(-log10(double(gmat{i})/65535)));
        b(i) = mean(mean(-log10(double(bmat{i})/65535)));
        x(i) = mean(mean(xmat{i}));
        y(i) = mean(mean(ymat{i}));
        z(i) = mean(mean(zmat{i}));
    end
        %%
    % figure;
    % plot(r,x,r,y,r,z);
    %%
    for i= 1 :nbbands
        X = xmat{i};
        X = X/mean(X(:))
        Z = zmat{i};
        T = xmat{i}./zmat{i};
        T = T/mean(T(:));
    %     %do not filter, it doesn't look good
    %     N = 1;
    %     f = ones(N,N)/N^2;
    %     figure;
    %     imagesc([conv2(X,f,'valid') conv2(T,f,'valid')]);
    %     colormap('gray');
    %     impixelinfo
        figure;
        imagesc([X T]);
        colormap('gray');
        impixelinfo
    end

    %%
    X = []; Y = []; Z = []; 
    for i= 1 :nbbands
       if strcmp(stackdir,'Vertical')
            k = 1;
       elseif strcmp(stackdir,'Horizontal')
            k = 2;   
       end
       X = cat(k,X,xmat{i});
       Y = cat(k,Y,ymat{i});
       Z = cat(k,Z,zmat{i});                           
    end
    %HB 11 July 2022: this was missing
    %the range is important to assure that the ECR do no go being 4.81
    %which corresponds to a signal of 1
    %znorm = mean(Z(:));
    %Z = Z/znorm;
    %HB 12 July: I want to try subtracting THETA by a min(trange) instead to adjust the range of THETA values.
    %
    figure;
    subplot(1,3,1);
    imagesc(X);
    set(gca,'DataAspectRatio',[1 1 1]);
    impixelinfo;
    colormap('gray'); colorbar;
    subplot(1,3,2);
    imagesc(Z);
    set(gca,'DataAspectRatio',[1 1 1]);
    impixelinfo;
    colormap('gray'); colorbar;
    subplot(1,3,3);
    imagesc(X./Z);
    set(gca,'DataAspectRatio',[1 1 1]);
    impixelinfo;
    colormap('gray'); colorbar;
    %signal correction range
    trange = [min(X(:)./Z(:)) max(X(:)./Z(:))];
    %% Writing the correcion file
    [ofilename,opathname] = uiputfile({'*.3ch'},'.lin signal calibration');
    if ofilename==0            
    else
        filename = fct_makecleanfilename(opathname,ofilename);
%         err = fct_WriteMultiCorrection(filename,Rot,rrange,grange,brange,trange,znorm);
%         [Rotcheck,rrange_,grange_,brange_,trange,znorm,err_]  = fct_ReadMultiCorrection(filename);
        err = fct_WriteMultiCorrection(filename,Rot,rrange,grange,brange,trange);
        [Rotcheck,rrange_,grange_,brange_,trange,err_]  = fct_ReadMultiCorrection(filename);
        Rotcheck-Rot
        E = Rot';
        Echeck = Rotcheck';
        %
        figure;
        hold on;
        plot3([0 1],[0 0],[0 0],'r','linewidth',2);
        plot3([0 0],[0 1],[0 0],'g','linewidth',2);
        plot3([0 0],[0 0],[0 1],'b','linewidth',2);
        plot3([-E(1,1) E(1,1)]/2,[-E(2,1) E(2,1)]/2,[-E(3,1) E(3,1)]/2,':','color',[1 1 0],'linewidth',2);
        plot3([-E(1,2) E(1,2)]/2,[-E(2,2) E(2,2)]/2,[-E(3,2) E(3,2)]/2,':','color',[0 1 1],'linewidth',2);
        plot3([-E(1,3) E(1,3)]/2,[-E(2,3) E(2,3)]/2,[-E(3,3) E(3,3)]/2,':','color',[1 0 1],'linewidth',2);
        plot3([-Echeck(1,1) Echeck(1,1)]/2,[-Echeck(2,1) Echeck(2,1)]/2,[-Echeck(3,1) Echeck(3,1)]/2,'--','color',[1 1 0],'linewidth',2);
        plot3([-Echeck(1,2) Echeck(1,2)]/2,[-Echeck(2,2) Echeck(2,2)]/2,[-Echeck(3,2) Echeck(3,2)]/2,'--','color',[0 1 1],'linewidth',2);
        plot3([-Echeck(1,3) Echeck(1,3)]/2,[-Echeck(2,3) Echeck(2,3)]/2,[-Echeck(3,3) Echeck(3,3)]/2,'--','color',[1 0 1],'linewidth',2);    
        hold off;
        legend('Red','Green','Blue','E#1','E#2','E#3');
        xlabel('x'); ylabel('y'); zlabel('z');    
    end
end
