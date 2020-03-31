function err = fct_CreateCalCurveSingleNewUncertROI(handles);
% 
% clear all;
% clear functions;
% close all;
% clc;
% fct_AddGafguiFctPath();
% handles = fct_initGafgui();

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
    [CAL_BCK,RES_BCK] = fct_read_tif_image(fct_makecleanfilename(ipathname1,ifilename1),'All');
    figure('NumberTitle','off','Name','Calibration curve: background');
    if length(size(CAL_BCK))==3
        RGB = 1;
    else
        RGB = 0;
    end
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
        if RGB==1
            Ibck = -log10(max(double(CAL_BCK(:,:,1)),1)./(2^(16)-1));
            Ibck = cat(3,Ibck,-log10(max(double(CAL_BCK(:,:,2)),1)./(2^(16)-1)));
            Ibck = cat(3,Ibck,-log10(max(double(CAL_BCK(:,:,3)),1)./(2^(16)-1)));
        else
            Ibck = -log10(max(double(CAL_BCK(:,:,1)),1)./(2^(16)-1));
        end

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
            if RGB==1
                BCKGRND = zeros(nbfilms,3);
            else
                BCKGRND = zeros(nbfilms,1);
            end

            NBCK = width_bck*height_bck/(handles.CCDres^2);

            figure(h);
            for i=1:nbfilms
                [newz,center,owidth,cropz] = fct_getroi(Ibck(:,:,1),RES_BCK,type,[width_bck height_bck]);
                [mu, sigmamu, s] = fct_analyze_region(newz);
                BCKGRND(i,1) = mu;
                if RGB==1
                    newz = imcrop(Ibck(:,:,2),cropz);
                    [mu, sigmamu, s] = fct_analyze_region(newz);
                    BCKGRND(i,2) = mu;
                    newz = imcrop(Ibck(:,:,3),cropz);
                    [mu, sigmamu, s] = fct_analyze_region(newz);
                    BCKGRND(i,3) = mu;
                end
            end
            figure(h);
            close(h);
            depthbck = length(size(CAL_BCK));
            clear CAL_BCK;
            clear Ibck;
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            [CAL_RAD,RES_RAD] = fct_read_tif_image(fct_makecleanfilename(ipathname2,ifilename2),'All');
            if length(size(CAL_RAD))~=depthbck
                error('Mismatch in image depths.');
            else
                if length(size(CAL_RAD))==3
                    RGB = 1;
                else
                    RGB = 0;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%
                figure('NumberTitle','off','Name','Courbe de calibration: films irradies');
                h = fct_display(CAL_RAD,RES_RAD);
                figure(h);
                ans = questdlg('Do you want to zoom to a smaller ROI?','Image','Yes','No','Yes') ;
                if strcmp(ans,'Yes')
                    figure(h);
                    CAL_RAD = imcrop;
                    h = fct_display(CAL_RAD,RES_RAD);
                end
                figure(h);

                if RGB==1
                    Irad = -log10(max(double(CAL_RAD(:,:,1)),1)./(2^(16)-1));
                    Irad = cat(3,Irad,-log10(max(double(CAL_RAD(:,:,2)),1)./(2^(16)-1)));
                    Irad = cat(3,Irad,-log10(max(double(CAL_RAD(:,:,3)),1)./(2^(16)-1)));
                else
                    Irad = -log10(max(double(CAL_RAD(:,:)),1)./(2^(16)-1));
                end

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
                        xrad(i) = center(1);
                        yrad(i) = center(2);
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
    else
        flag = 0;
    end
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
    %         button = questdlg('Which channel do you want to analyze?','Channel','Red','Green','Blue','Red');
    %     channel = fct_colortochannel(button)
    nchannel = size(Irad,3);
    for channel =1:nchannel;
        if nchannel>1
            button = questdlg(['Do you want to characterize the ' sprintf('%s',fct_channeltocolor(channel)) ' channel?'],'Mutichannel correction','Yes','No','Yes') ;
        else
            button = 'Yes';
        end
        if strcmp(button,'Yes')

            Iradtmp = Irad(:,:,channel);
            BCKGRNDtmp = BCKGRND(:,channel);
            
%             button = questdlg('What kind of background analysis do you want to do?','Background analysis','Average','One-by-one','Average') ;
            button = 'Average';
            if strcmp(button,'Average')
                %ici ajouter un while au cas ou cela ne marche pas
                BCKGRNDtmp = BCKGRNDtmp*0+mean(BCKGRNDtmp);
            end
           
            
            DOSE = DOSE(:);
            [DOSE,isort] = sort(DOSE);
            RAWOD = zeros(nbfilms,1);

            [DOSE,isort] = sort(DOSE);
            numpix = [];
            A = [];
            XI = []; sXI = [];
            xi = []; muxi = [];
            dose = [];
            for i=isort(:)' %this make them in order of dose instead of i=1:nbfilms
                [nlines,ncols] = size(Iradtmp);
                [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,RES_RAD);
                xindex = fct_postoindex(xrad(i)-width_rect/2,xgrid);
                yindex = fct_postoindex(yrad(i)-height_rect/2,ygrid);
                A{i} = imcrop(Iradtmp,[xindex yindex width_rect/RES_RAD height_rect/RES_RAD]);
                [mu, sigma, s] = fct_analyze_region(A{i});
                mu = mu-BCKGRNDtmp(i);
                [m,n] = size(A{i});
                XI = cat(1,XI,mu);
                sXI = cat(1,sXI,sigma);
                xi = cat(1,xi,reshape(A{i},m*n,1));
                muxi = cat(1,muxi,reshape(A{i},m*n,1)*0+mu);
                dose = cat(1,dose,reshape(A{i},m*n,1)*0+DOSE(i));
                numpix = cat(1,numpix,n*m);
            end
            OD = XI(isort);
            %%%backup for autoconsitency test
            ODinit = OD;
            DOSEinit = DOSE;
            Npix = width_rect*height_rect/handles.CCDres/handles.CCDres;
            
            j=0;
            k = 0;
            option = [1 2 4];
            Nrange = 2:min((nbfilms-2),10);
            
            hwait = waitbar(1/length(Nrange)/length(option),'Best fit calculation');
            
            optval = [];
            Nval = [];
            h_fig = [];
            sycal = [];
            Errdose = [];
            
            for opt = option
                SY = [];
                for N = Nrange
                    [p,sy,R,df] = fct_lsf(DOSE,OD,N,opt);
                    if fct_iscurvevalid(p,N,opt,min(DOSE),max(DOSE))
                        SY = [SY sy];
                    else
                        SY = [SY 10^6];
                    end
                    j = j+1;
                    waitbar(j/length(Nrange)/length(option),hwait);
                    figure(hwait);
                end
                sy = min(SY);
                ii = min(find(SY==sy));
                N = Nrange(ii);
                F = fct_F_matrix(DOSE,N,opt);
                [p,sy,R,df] = fct_lsf(DOSE,OD,N,opt);
                y = F*p;
                test = norm(y - sort(y));
                %[k,sy,N,opt]
                [htmp,stmp,errtmp] = fct_show_calcurve(DOSE,OD,N,opt,0);
                if htmp~=-1
                    k = k+1;
                    h_fig{k} = htmp;
                    sycal{k} = stmp;
                    Errdose{k} = errtmp;                    
                    Nval{k} = N;
                    optval{k} = opt;
                end                    
            end            
            close(hwait);
            
%             %This code no longer works with 2014b (perhaps even before): #£"$%$&!!
%             h_msg = helpdlg('Choose curve and press enter','Calibration curve');
%             uiwait(h_msg) ;                        
%             w = 0;
%             figure(handles.H);
%             while(w==0)
%                 w = waitforbuttonpress
%             end

            s = zeros(k,1);
            for i=1:k
                s(i) = sycal{i};
            end
            k = min(find(s==min(s)));
            
            h = figure(h_fig{k});
            for i = 1:k
                if h == h_fig{i}
                    K = i;
                else
                    if ishandle(h_fig{i})
                        close(h_fig{i});
                    end
                end
            end
            h = h_fig{K};
            N = Nval{K};
            opt = optval{K};
            [p,sycal,R,df] = fct_lsf(DOSE,OD,N,opt);
            
            Nsave = N;
            optsave = opt;
            Fsave = fct_F_matrix(DOSE,Nsave,optsave);
            [psave,sy,R,df] = fct_lsf(DOSE,OD,Nsave,optsave);
            XIhat = Fsave*psave;             
%             figure;
%             plot(DOSE,OD,'or',DOSE,XIhat,'xk');            
            npix0 = (RES_RAD)^2/handles.CCDres^2;
            %%%%%%
            %Here we need to evaluate how the uncertainty behaves with ROI    
            npix = [];    
            sxi = []; 
            sXI = [];
            NB = [];
            for f = 2.^[0:4];
                a = fct_reducematrix(A{1},f);
                [m,n] = size(a);
                if (n*m)>20
                    sxi_reduce = zeros(length(DOSE),1);
                    npix = cat(1,npix,(RES_RAD*f)^2/handles.CCDres^2);                    
                    nb = zeros(length(DOSE),1);
%                     figure; hold on;
                    for k = 1:length(DOSE)
                        a = fct_reducematrix(A{k},f);
                        [m1,n1] = size(a);
                        xi_reduce = reshape(a,m1*n1,1)-BCKGRNDtmp(k);
                        sxi_reduce(k) = sqrt(sum((xi_reduce-XIhat(k)).^2)/(m1*n1-Nsave));
                        nb(k) = n1*m1;
%                         plot(DOSE(k)+0*xi_reduce,xi_reduce,'.k'); 
                    end
%                     hold off;
%                     sXI = cat(2,sXI,sqrt(sum(sxi_reduce.^2.*(nb-1))/(sum(nb)-Nsave)));
                    sxi = cat(2,sxi,sxi_reduce(:));
                    sXI = cat(1,sXI,sqrt(sum(sxi_reduce.^2.*(nb-1))/(sum(nb)-Nsave)));
                    NB = cat(2,NB,nb(:));
                end
            end
            
            figure; 
            subplot(1,2,1);
            hold on;
            for i=1:length(npix)
                plot(DOSE,sxi(:,i));
                tag{i} = sprintf('%.2f x %.2f mm^2',10*sqrt(npix(i))*handles.CCDres,10*sqrt(npix(i))*handles.CCDres);
            end
            hold off;
            xlabel('Dose');
            ylabel('sOD');
            legend(tag,'location','best');
%             y = ylim;
%             ylim([0 y(2)]);
            subplot(1,2,2);
            hold on;
            x = 1./npix;
            y = sXI.^2;
            p = polyfit(x,y,1);
            plot(1./npix,sXI.^2,'.k',x,polyval(p,x),'--r');
            hold off;
            xlabel('1/Npix');
            ylabel('sOD^2');
            
            sig0 = sqrt(p(2));
            sig1 = sqrt(p(1));
            
%             else
%                 sig0 = sycal/sqrt(2);
%                 sig1 = 0;
%             end
            
            [ofilename,opathname]=uiputfile({'*.cal'},'Save calibration curve');
            size(DOSE)
            size(OD)
            
            if ~strcmp(class(ofilename),'double')
                file = fopen(fct_makecleanfilename(opathname,ofilename),'w');
                for i=1:nbfilms
                    fprintf(file,'%f\t%f\n',DOSE(i),OD(i));
                end
                %HB 15 April 2015: this is the new way of dealing with background
                %we take the average during calibration and we assume the ROI will be
                %large enough during use        
                NBCK = 1e10;
                if strcmp(button,'Average')
                    sbck = std(BCKGRND(:,channel));
                    sig0 = sqrt(sig0.^2+sbck.^2)/sqrt(2);
                else
                    sig0 = sig0/sqrt(2);
                end
                fprintf(file,'%d\t%d\n',Npix,NBCK);
                fprintf(file,'%f\t%f\n',sig0,sig1);
                %fprintf(file,'%f\t%f\n',sycal,(height_rect/handles.CCDres)*(width_rect/handles.CCDres));
                fprintf(file,'%d\t%d',optsave,Nsave);
                fclose(file);
            end
            close(h);
        end
    end
end

err = 1;