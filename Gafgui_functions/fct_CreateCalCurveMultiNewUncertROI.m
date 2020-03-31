% function err = fct_CreateCalCurveMultiNewUncertROI(handles);

clear all;
clear functions;
close all;
clc;
fct_AddGafguiFctPath();
handles = fct_initGafgui('Black');

err = 0;

button = questdlg('Do you want to start from the beginning?','Shortcut','Yes','No','Yes') ;
handles = fct_updatedisplay(handles);
figure(handles.H);

if strcmp(button,'No')
    flag = 1;
else
    [ifilename1,ipathname1]=uigetfile({'*.tif'},'Choose unirradiated films scan');
    if (~strcmp(class(ifilename1),'double'))
        [ifilename2,ipathname2]=uigetfile({'*.tif'},'Choose irradiated films scan');
        if (~strcmp(class(ifilename2),'double'))
            flag = 2;
        else
            flag = 0;
        end
    else
        flag = 0;
    end
end

if flag == 2      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [CAL_BCK,RES_BCK,BITS_BCK,CHANNEL_BCK] = fct_read_tif_image(fct_makecleanfilename(ipathname1,ifilename1),'All');
    [CAL_RAD,RES_RAD,BITS,CHANNEL] = fct_read_tif_image(fct_makecleanfilename(ipathname2,ifilename2),'All');
    
    if (BITS_BCK~=16)||(CHANNEL_BCK~=4)||(BITS~=16)||(CHANNEL~=4)
        errordlg('Wrong type of tiff. Try grayscale 16 bits.');
        flag = 0;
    else
        figure('NumberTitle','off','Name','Dose calibration: unirradiated films');
        h = fct_display(CAL_BCK,RES_BCK);

        ans = questdlg('Do you want to zoom to a smaller ROI?','Image','Yes','No','Yes') ;
        if strcmp(ans,'Yes')
            figure(h);
            CAL_BCK = imcrop;
            h = fct_display(CAL_BCK,RES_BCK);
        end
        figure(h);
        nbfilms  = inputdlg({'Number of films to calibrate'},'Number of films',1);
        if numel(nbfilms)~=0
            nbfilms  = str2double(nbfilms);

            Ibck = -log10(max(double(CAL_BCK(:,:,1)),1)./(2^(16)-1));

            str{1}='Fixed';
            str{2}='Zoom-fixed';
            str{3}='Point';

            [type,ok] = listdlg('Name','ROI select type','ListString',str);

            if ok==1
                type = str{type};
            else
                type = 'Fixed';
            end

            rectsize  = inputdlg({'Region width in cm:','Region height in cm:'},'Background',1,{'1','1'}) ;
            if numel(rectsize)==0
                flag = 0;
            else
                width_bck = str2double(rectsize(1));
                height_bck = str2double(rectsize(2));
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                BCKGRND = zeros(nbfilms,1);

                NBCK = width_bck*height_bck/(handles.CCDres^2);

                figure(h);
                for i=1:nbfilms
                    [newz,center,owidth,cropz] = fct_getroi(Ibck(:,:,1),RES_BCK,type,[width_bck height_bck]);
                    [mu, sigmamu, s] = fct_analyze_region(newz);
                    BCKGRND(i,1) = mu;
                end
                figure(h);
                close(h);
                clear CAL_BCK;
                clear Ibck;
                %%%%%%%%%%%%%%%%%%%%%%%%%%

                figure('NumberTitle','off','Name','Dose calibration: irradiated films');
                h = fct_display(CAL_RAD,RES_RAD);
                figure(h);
                ans = questdlg('Do you want to zoom to a smaller ROI?','Image','Yes','No','Yes') ;
                if strcmp(ans,'Yes')
                    figure(h);
                    CAL_RAD = imcrop;
                    h = fct_display(CAL_RAD,RES_RAD);
                end
                figure(h);

                Irad = -log10(max(double(CAL_RAD(:,:,1)),1)./(2^(16)-1));

                rectsize  = inputdlg({'Region width in cm:','Region height in cm:'},'Irradiated films',1,{'1','1'}) ;
                if numel(rectsize)==0
                    flag = 0;
                else
                    width_rect = str2double(rectsize(1));
                    height_rect = str2double(rectsize(2));

                    [type,ok] = listdlg('Name','ROI select type','ListString',str);

                    if ok==1
                        type = str{type};
                    else
                        type = 'Fixed';
                    end
                    figure(h);
                    xrad(1:nbfilms) = 0;
                    yrad(1:nbfilms) = 0;
                    for i=1:nbfilms
                        [newz,center,owidth] = fct_getroi(Irad(:,:,1),RES_RAD,type,[width_rect height_rect]);
                        xrad(i) = center(1)
                        yrad(i) = center(2)
                    end
                    figure(h);
                    close(h);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%

                    DOSE = [];
                    nprompt = ceil(nbfilms/5);

                    button = questdlg('How do you wish to enter dose values?','Dose values','File','Manual','File') ;

                    flag = 1;
                    if strcmp(button,'Manual')
                        for j = 1:nprompt
                            clear prompt;
                            for i= 1+5*(j-1):min(5*j,nbfilms)
                                dumb = sprintf('Film #%d\nDose (CMU)',i);
                                prompt{i-5*(j-1)} = dumb;
                            end
                            answer = inputdlg(prompt,'Dose values (CMU)',1);
                            DOSE = cat(1,DOSE,str2double(answer));
                        end
                        DOSE = DOSE';
                        %ici il devrait y avoir une facon de valider les valeurs de dose
                        DOSE = abs(DOSE);
                    else
                        [ifilename,ipathname] = uigetfile({'*.txt'},'Choose file containing dose values');
                        if ~strcmp(class(ifilename),'double')
                            file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');
                            DOSE = fscanf(file,'%f',[1 inf]);
                        else
                            flag = 0;
                        end
                    end
                    if flag~=0
                        button = questdlg('Do you want to save this step?','Shortcut','Yes','No','Yes') ;

                        if strcmp(button,'Yes')
                            %ici ajouter un while au cas ou cela ne marche pas
                            [ofilename,opathname]=uiputfile({'*.mat'},'Work to save');
                            filename = fct_makecleanfilename(opathname,ofilename);
                            save(filename,'Irad','xrad','yrad','nbfilms','RES_RAD','NBCK','BCKGRND','DOSE','width_rect','height_rect');
                        end

                        clear CAL_RAD;

                        figure(handles.H);
                    end
                end
            end
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flag == 1
    
    [ifilename,ipathname]=uigetfile({'*.mat'},'Work to import');
    if ~strcmp(class(ifilename),'double')        
        filename = fct_makecleanfilename(ipathname,ifilename);
        load (filename);
    else
        flag = 0;
    end
end

if flag~=0

    channel = 1;
    Iradtmp = Irad(:,:,channel);

    XI = [];
    sXI = [];
    xi = [];
    muxi = [];
    dose = [];

    [DOSE,isort] = sort(DOSE);
    numpix = [];
    for i=isort %this make them in order of dose instead of i=1:nbfilms
        [nlines,ncols] = size(Iradtmp);
        [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,RES_RAD);
        xindex = fct_postoindex(xrad(i)-width_rect/2,xgrid);
        yindex = fct_postoindex(yrad(i)-height_rect/2,ygrid);
        A{i} = imcrop(Iradtmp,[xindex yindex width_rect/RES_RAD height_rect/RES_RAD]);
        [mu, sigma, s] = fct_analyze_region(A{i});
        [m,n] = size(A{i});
        XI = cat(1,XI,mu);
        sXI = cat(1,sXI,sigma);
        xi = cat(1,xi,reshape(A{i},m*n,1));
        muxi = cat(1,muxi,reshape(A{i},m*n,1)*0+mu);
        dose = cat(1,dose,reshape(A{i},m*n,1)*0+DOSE(i));
        numpix = cat(1,numpix,n*m);
    end

    G = [DOSE(:).^0./sXI(:) DOSE(:)./sXI(:)];
    z = [XI(:)./sXI(:)];
    c = inv(G'*G)*G'*z;
    XIhat = [DOSE(:).^0 DOSE(:)]*c;
    
    sxi = [];
    for j=1:length(DOSE)
        k = find(dose==DOSE(j));
        tmp = sqrt(sum((xi(k)-XIhat(j)).^2)/(length(k)-2));
        sxi = cat(1,sxi,tmp);
    end
    sXI = sxi;
    
    %account for background value
    %this is temporary to fix a bug
    BCKGRND = BCKGRND(:);
    k = find(BCKGRND~=0);
    BCKGRND = BCKGRND(k);
    
    sbck = std(BCKGRND);
    
    %Uncertainty vs Npix
    %The default resolution first
    npix0 = RES_RAD^2/handles.CCDres^2;
%     F = [XI(:).^0 (XI(:)-0).^1 (XI(:)-0).^2];
    F = [XI(:).^0 (XI(:)-0.1).^2];
    coeff = inv(F'*F)*F'*sXI(:);
    sXIhat0 = F*coeff;
    W = diag(sXIhat0.^2);
    V = (inv(G'*G)*G')*W*(inv(G'*G)*G')';
    figure;
    plot(XI,sXI,'xm',XI,sXIhat0,'k','linewidth',2);
    title('Uncertainty on xi');
    xlabel('XI'); ylabel('sXI');  exp = ['^2'];
    str{1} = sprintf('Observed: res = (%.2f mm)%s',RES_RAD*10,exp);
    str{2} = sprintf('Fitted: res = (%.2f mm)%s',RES_RAD*10,exp);
    legend(str,'location','northwest');
    
    %Here we need to evaluate how the uncertainty behaves with ROI    
    npix = [];    
    sxi = [];
    val = [];
    figure; hold on;
    for f = 2.^[0:4];
        a = fct_reducematrix(A{1},f);
        [m,n] = size(a);
        if (n*m)>50
            sxi_reduce = zeros(length(DOSE),1);
            npix = cat(1,npix,(RES_RAD*f)^2/handles.CCDres^2);
            for k = 1:length(DOSE)
                [m1,n1] = size(A{k});
                b = reshape(A{k},m1*n1,1);
                %Italian statistics =)
%                 [~,i] = sort(rand(m1*n1,1));
%                 i = 1:(m1*n1);
%                 b = reshape(b(i),m1,n1);
                a = fct_reducematrix(A{k},f);
                [m1,n1] = size(a);
                xi_reduce = reshape(a,m1*n1,1);
                tmp = sqrt(sum((xi_reduce-XIhat(k)).^2)/(length(xi_reduce)-2));
%                 tmp = std(xi_reduce);
                sxi_reduce(k) = tmp;
            end
            coeff = inv(F'*F)*F'*sxi_reduce(:);
            sxi = cat(2,sxi,F*coeff);     
%             plot(XI,sxi_reduce,'.r',XI,(F*coeff,'k');
        end
    end
%     hold off;
%     xlabel('Xi');
%     ylabel('sXi');

    %seperate measurements
    bck =[  0.0739 0.0068
            0.0753 0.0082
            0.0694 0.0068
            0.076 0.0069
            0.0712 0.0077
            0.0722 0.0076];
    sbck = sqrt(std(bck(:,1))^2+max(bck(:,2))^2);    
    [m,n] = size(sxi);
%     SXI = sxi/sxi(1,n)*sqrt(sxi(1,n)^2+sbck^2); 
    SXI = sxi/sxi(1,1)*sqrt(sxi(1,1)^2+sbck^2); 
    sxi = SXI;

    [ni,xi] = meshgrid(npix(:),XI(:));
    figure;
    mesh(sqrt(npix0./ni),xi,sxi);
    xlabel('sqrt(N_0/N_1)');
    ylabel('xi');
    zlabel('sxi');
    
        
    sximesh.x = sqrt(npix0./ni);
    sximesh.y = xi;
    sximesh.z = sxi;
    sximesh.npix0 = npix0;
    %Finally, estimate for 1 cm^2
    npix = floor(0.1/handles.CCDres^2);
    sXIhat = fct_GetsXI(sximesh,npix,XI);
    %Now plot it
    figure;
    subplot(1,2,1);
    plot(XI,sXI,'xm',XI,sXIhat0,'k',XI,sXIhat,'b','linewidth',2);
    title('Uncertainty on xi versus dose');
    xlabel('XI'); ylabel('sXI');
    exp = ['^2'];
    str{1} = sprintf('Observed: res = (%.2f mm)%s',RES_RAD*10,exp);
    str{2} = sprintf('Fitted: res = (%.2f mm)%s',RES_RAD*10,exp);
    str{3} = sprintf('Predicted: : res = 1 mm%s',exp);
    legend(str);
    
    %Now uncertainty on dose for a given ROI
    dose = [min(DOSE):2:max(DOSE)]';
    npix = floor(0.1/handles.CCDres^2);
    sd1 = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,npix0,dose);
    sd2 = fct_UncertaintyMultichannelSingleMeas(c,V,sximesh,npix,dose);    
    subplot(1,2,2);
    plot(dose,sd1./dose(:)*100,'k',dose,sd2./dose(:)*100,'b','linewidth',2);
    ylabel('Dose uncertainty (%)');
    ylim([0 5]);
    title('Uncertainty on dose for single measurements');
    xlabel('Dose index');
    exp = ['^2'];
    str{1} = sprintf('Predicted: res = (%.2f mm)%s',RES_RAD*10,exp);
    str{2} = sprintf('Predicted: res = 1 mm%s',exp);
    legend(str);
    
    h = fct_show_calcurve_multi(DOSE,XI,sXI,c,V,RES_RAD,sximesh);    
    [m,n] = size(ni);
    
    [ofilename,opathname]=uiputfile({'*.mlt'},'Save calibration curve');
    if ~strcmp(class(ofilename),'double')
        file = fopen(fct_makecleanfilename(opathname,ofilename),'w');
        fprintf(file,'%.10e\t%.10e\t%.10e\n',c(1),c(2),RES_RAD);
        fprintf(file,'%.10e\t%.10e\t%.10e\n',V(1,1),V(1,2),V(2,2));       
        fprintf(file,'%.10e\t%.10e\t%.10e\n',m,n,sximesh.npix0);
        for j=1:n
            for i=1:m
                fprintf(file,'%.10e\t%.10e\t%.10e\n',sximesh.x(i,j),sximesh.y(i,j),sximesh.z(i,j));                
            end
        end
        for i=1:length(DOSE)-1
            fprintf(file,'%.10e\t%.10e\t%.10e\n',DOSE(i),XI(i),sXIhat0(i));
        end
        fprintf(file,'%.10e\t%.10e\t%.10e',DOSE(i+1),XI(i+1),sXIhat0(i+1));
        fclose(file);
    end    
end

err = 1;