function handles = fct_CreateSignalLinearizationMap(handles)
%% 
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
    end
end
if flag==2%data acquisition of dark and bright films
    %%%
    [ifilename,ipathname]=uigetfile({'*.tif'},'Choose the dark scan');
    if ~strcmp(class(ifilename),'double')
        [Iout,RES] = fct_read_tif_image(fct_makecleanfilename(ipathname,ifilename),'All');
        smin(1) = floor(mean(mean(Iout(:,:,1))));
        smin(2) = floor(mean(mean(Iout(:,:,2))));
        smin(3) = floor(mean(mean(Iout(:,:,3))));
    else
        htmp = msgbox('No dark scan selected, minimal signal set to 0.');
        uiwait(htmp);
        smin = [0 0 0];
        %flag = -1;
    end
    [ifilename,ipathname]=uigetfile({'*.tif'},'Choose the bright scan');
    if ~strcmp(class(ifilename),'double')
        [Iout,RES] = fct_read_tif_image(fct_makecleanfilename(ipathname,ifilename),'All');
        smax(1) = ceil(mean(mean(Iout(:,:,1))));
        smax(2) = ceil(mean(mean(Iout(:,:,2))));
        smax(3) = ceil(mean(mean(Iout(:,:,3))));
    else
        htmp = msgbox('No dark scan selected, maximal signal set to 65535.');
        uiwait(htmp);
        smin = [65535 65535 65535];
        %flag = -1;
    end
    smin = [min(smin(1),smax(1)) min(smin(2),smax(2)) min(smin(3),smax(3))]
    smax = [max(smin(1),smax(1)) max(smin(2),smax(2)) max(smin(3),smax(3))]
end
if flag==2%data acquisition of NIST OD values
    %%%
    [ifilename,ipathname] = uigetfile({'*.txt'},'Choose file with NIST OD values');
    if ~strcmp(class(ifilename),'double')
        file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');
    else
        htmp = msgbox('NIST OD values not set. File not read.');
        uiwait(htmp);
        flag = -1;
    end
    if flag==2
        A = fscanf(file,'%f',[1 inf]);
        A = A';
        ODcal = A(:);
        Tcal = 10.^(-ODcal);
        nfilms = numel(rmat);
        if length(Tcal)<nfilms
            error('NIST OD values not set. Not enough NIST OD values for the number of films pieces.');
            flag = -1;
        end
        if flag==2
            Tcal = Tcal(1:nfilms);
            Tcal = sort(Tcal);

            mux = zeros(nfilms,1);muy = zeros(nfilms,1);muz = zeros(nfilms,1);
            sx = zeros(nfilms,1); sy = zeros(nfilms,1); sz = zeros(nfilms,1);
            for j=1:nfilms
                x = double(rmat{j}); 
                y = double(gmat{j}); 
                z = double(bmat{j}); 
                mux(j) = mean(x(:)); muy(j) = mean(y(:)); muz(j) = mean(z(:));
                sx(j) = std(x(:)); sy(j) = std(y(:)); sz(j) = std(z(:));
            end
            [~,j] = sort(mux);
            mux = mux(j); sx = sx(j);
            muy = muy(j); sy = sy(j);
            muz = muz(j); sz = sz(j);
        end
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
%                 save(filename,'posx','posy','rmat','bmat','gmat','bits','res');
            save(filename,'mux','muy','muz','Tcal','smin','smax');
        end
    end
end

%% Calibration
if flag ~=-1
    
    button = questdlg('Do you want a black-and-white only LUT?','B&W','Yes','No','No') ;
    if strcmp(button,'Yes')
        smin = ones(3,1)*mean(smin);
        smax = ones(3,1)*mean(smax);
        mux = mean([mux(:) muy(:) muz(:)],2);
        muy = mux;
        muz = muy;
    end
    
    calorder = 4;

    x = (mux-smin(1))/(smax(1)-smin(1));
    y = (muy-smin(2))/(smax(2)-smin(2));
    z = (muz-smin(3))/(smax(3)-smin(3));
    Fx = []; Fy = []; Fz = []; 
    %
    for n=1:calorder
            Fx = cat(2,Fx,x(:).^n);
            Fy = cat(2,Fy,y(:).^n);
            Fz = cat(2,Fz,z(:).^n);
    end
    %fits
    pcal{1} = Fx\Tcal;
    pcal{2} = Fy\Tcal;
    pcal{3} = Fz\Tcal;
    %estimations
    xhat = Fx*pcal{1};
    yhat = Fy*pcal{2};
    zhat = Fz*pcal{3};
    % %This makes pcal usable with polyval
    % pcal{1} = [flipud(pcal{1}(:))' 0];
    % pcal{2} = [flipud(pcal{2}(:))' 0];
    % pcal{3} = [flipud(pcal{3}(:))' 0];
    % xhat = polyval(pcal{1},(mux-smin(1))/smax(1));
    % yhat = polyval(pcal{2},(muy-smin(2))/smax(2));
    % zhat = polyval(pcal{3},(muz-smin(3))/smax(3));
    %
    figure;
    subplot(1,2,1);
    hold on;
    plot(mux,Tcal,'or','linewidth',2); 
    plot(muy,Tcal,'og','linewidth',2); 
    plot(muz,Tcal,'ob','linewidth',2);
    plot(mux,xhat,'--r','linewidth',2); 
    plot(muy,yhat,'--g','linewidth',2); 
    plot(muz,zhat,'--b','linewidth',2);
    legend('Red channel','Green channel','Blue channel','location','northwest');
    hold off;
    xlabel('Raw signal');
    ylabel('Transmittance');

    subplot(1,2,2);
    hold on;
    plot(mux,xhat-Tcal,'--or','linewidth',2); 
    plot(muy,yhat-Tcal,'--og','linewidth',2); 
    plot(muz,zhat-Tcal,'--ob','linewidth',2);  
    legend('Red channel','Green channel','Blue channel','location','northwest');
    hold off;
    xlabel('Raw signal');
    ylabel('Residual of transmittance');
    %% Make Look Up Table (LUT)
    %max values of normalized r,g,b
    irmax = (max(mux)-smin(1))/(smax(1)-smin(1));
    igmax = (max(muy)-smin(2))/(smax(2)-smin(2));
    ibmax = (max(muz)-smin(3))/(smax(3)-smin(3));
    %map of normalized r,g,b for LUT
    dr = 1;%(smax(1)-smin(1))/(2^16-1);
    dg = 1;%(smax(2)-smin(2))/(2^16-1);
    db = 1;%(smax(3)-smin(3))/(2^16-1);
    ir = ([smin(1):dr:smax(1)]'-smin(1))/(smax(1)-smin(1));
    ig = ([smin(2):dg:smax(2)]'-smin(2))/(smax(2)-smin(2));
    ib = ([smin(3):db:smax(3)]'-smin(3))/(smax(3)-smin(3));
    Fir = []; Fig = []; Fib = []; 
    for n=1:calorder
            Fir = cat(2,Fir,ir(:).^n);
            Fig = cat(2,Fig,ig(:).^n);
            Fib = cat(2,Fib,ib(:).^n);
    end
    itr = Fir*pcal{1};
    itg = Fig*pcal{2};
    itb = Fib*pcal{3};
    % %end portion overwritten by linear
    % k = min(find(itr>=irmax));
    % itr(k:end) = interp1([ir(k) 1],[itr(k) 1],ir(k:end));
    % k = min(find(itg>=igmax));
    % itg(k:end) = interp1([ig(k) 1],[itg(k) 1],ig(k:end));
    % k = min(find(itb>=ibmax));
    % itb(k:end) = interp1([ib(k) 1],[itb(k) 1],ib(k:end));

    r = (smin(1) + ir*(smax(1)-smin(1)));
    g = (smin(2) + ig*(smax(2)-smin(2)));
    b = (smin(3) + ib*(smax(3)-smin(3)));
    rcorr = itr*(2^16-1);
    gcorr = itg*(2^16-1);
    bcorr = itb*(2^16-1);

    %extrapolate/clip
    rcorr = [zeros(1,length(0:(min(r)-1))) rcorr(:)' ones(1,length((max(r)+1):(2^16-1)))*max(rcorr)]';
    r     = [               0:(min(r)-1)   r(:)'                    (max(r)+1):(2^16-1)]';
    gcorr = [zeros(1,length(0:(min(g)-1))) gcorr(:)' ones(1,length((max(g)+1):(2^16-1)))*max(gcorr)]';
    g     = [               0:(min(g)-1)   g(:)'                    (max(g)+1):(2^16-1)]';
    bcorr = [zeros(1,length(0:(min(b)-1))) bcorr(:)' ones(1,length((max(b)+1):(2^16-1)))*max(bcorr)]';
    b     = [               0:(min(b)-1)   b(:)'                    (max(b)+1):(2^16-1)]';

    figure;
    hold on;
    plot(uint16(r),uint16(rcorr),'.r','linewidth',2); 
    plot(uint16(g),uint16(gcorr),'.g','linewidth',2); 
    plot(uint16(b),uint16(bcorr),'.b','linewidth',2);
    plot(mux,Tcal*65535,'or','linewidth',2); 
    plot(muy,Tcal*65535,'og','linewidth',2); 
    plot(muz,Tcal*65535,'ob','linewidth',2);
    % plot(ir,itr,'.r','linewidth',2); 
    % plot(ig,itg,'.g','linewidth',2); 
    % plot(ib,itb,'.b','linewidth',2);
    % plot([irmax irmax],[0 1],'--r','linewidth',2); 
    % plot([irmax irmax],[0 1],'--g','linewidth',2); 
    % plot([irmax irmax],[0 1],'--b','linewidth',2); 
    legend('Red channel','Green channel','Blue channel','location','northwest');
    hold off;
    xlabel('Raw signal');
    % ylabel('Transmittance');
    ylabel('Corrected signal');
end
%% Save calibration
if flag~=-1
    [ofilename,opathname] = uiputfile({'*.lin'},'.lin signal calibration');
    if ofilename==0            
    else
        filename = fct_makecleanfilename(opathname,ofilename);
        err = fct_WriteSignalCorrection(filename,rcorr,gcorr,bcorr);
%         %check
%         LUT = fct_ReadSignalCorrection(filename);
%         figure;
%         hold on;
%         plot(LUT.raw,LUT.red,'--r','linewidth',1); 
%         plot(LUT.raw,LUT.green,'--g','linewidth',1); 
%         plot(LUT.raw,LUT.blue,'--b','linewidth',1);
%         plot(uint16(r),uint16(rcorr),'or','markersize',2); 
%         plot(uint16(g),uint16(gcorr),'og','markersize',2); 
%         plot(uint16(b),uint16(bcorr),'ob','markersize',2);
%         legend('Red channel - LUT','Green channel - LUT','Blue channel - LUT','location','northwest');
%         hold off;
%         xlabel('Raw signal');
%         ylabel('Corrected signal');
    end
end