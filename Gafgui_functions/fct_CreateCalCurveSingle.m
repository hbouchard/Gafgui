function err = fct_CreateCalCurveSingle(handles);

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
%     CAL_BCK(:,:,1) = fliplr(CAL_BCK(:,:,1));
%     CAL_BCK(:,:,2) = fliplr(CAL_BCK(:,:,2));
%     CAL_BCK(:,:,3) = fliplr(CAL_BCK(:,:,3));
    figure('NumberTitle','off','Name','Calibration curve: background');
    h = fct_display(CAL_BCK,RES_BCK);
%     ans = questdlg('Do you want to flip the image horizontally?','Image','Yes','No','No') ;
%     if strcmp(ans,'Yes')
%         CAL_BCK(:,:,1) = fliplr(CAL_BCK(:,:,1));
%         CAL_BCK(:,:,2) = fliplr(CAL_BCK(:,:,2));
%         CAL_BCK(:,:,3) = fliplr(CAL_BCK(:,:,3));
%         close(h);
%         figure('NumberTitle','off','Name','Calibration curve: background');
%         h = fct_display(CAL_BCK,RES_BCK);
%     end
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
        Ibck = cat(3,Ibck,-log10(max(double(CAL_BCK(:,:,2)),1)./(2^(16)-1)));
        Ibck = cat(3,Ibck,-log10(max(double(CAL_BCK(:,:,3)),1)./(2^(16)-1)));

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
            BCKGRND = zeros(nbfilms,3);

            NBCK = width_bck*height_bck/(handles.CCDres^2);

            figure(h);
            for i=1:nbfilms
                [newz,center,owidth,cropz] = fct_getroi(Ibck(:,:,1),RES_BCK,type,[width_bck height_bck]);
                [mu, sigmamu, s] = fct_analyze_region(newz);
                BCKGRND(i,1) = mu;
                newz = imcrop(Ibck(:,:,2),cropz);
                [mu, sigmamu, s] = fct_analyze_region(newz);
                BCKGRND(i,2) = mu;
                newz = imcrop(Ibck(:,:,3),cropz);
                [mu, sigmamu, s] = fct_analyze_region(newz);
                BCKGRND(i,3) = mu;
            end
            figure(h);
            close(h);
            clear CAL_BCK;
            clear Ibck;
            %%%%%%%%%%%%%%%%%%%%%%%%%%

            [CAL_RAD,RES_RAD] = fct_read_tif_image(fct_makecleanfilename(ipathname2,ifilename2),'All');
%             CAL_RAD(:,:,1) = fliplr(CAL_RAD(:,:,1));
%             CAL_RAD(:,:,2) = fliplr(CAL_RAD(:,:,2));
%             CAL_RAD(:,:,3) = fliplr(CAL_RAD(:,:,3));
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            figure('NumberTitle','off','Name','Courbe de calibration: films irradies');
            h = fct_display(CAL_RAD,RES_RAD);
            figure(h);
%             ans = questdlg('Do you want to flip the image horizontally?','Image','Yes','No','No') ;
%             if strcmp(ans,'Yes')
%                 CAL_RAD(:,:,1) = fliplr(CAL_RAD(:,:,1));
%                 CAL_RAD(:,:,2) = fliplr(CAL_RAD(:,:,2));
%                 CAL_RAD(:,:,3) = fliplr(CAL_RAD(:,:,3));
%                 close(h);
%                 figure('NumberTitle','off','Name','Calibration curve: irradiated film');
%                 h = fct_display(CAL_RAD,RES_RAD);
%                 figure(h);
%             end
            ans = questdlg('Do you want to zoom to a smaller ROI?','Image','Yes','No','Yes') ;
            if strcmp(ans,'Yes')
                figure(h);
                CAL_RAD = imcrop;
                h = fct_display(CAL_RAD,RES_RAD);
            end
            figure(h);

            Irad = -log10(max(double(CAL_RAD(:,:,1)),1)./(2^(16)-1));
            Irad = cat(3,Irad,-log10(max(double(CAL_RAD(:,:,2)),1)./(2^(16)-1)));
            Irad = cat(3,Irad,-log10(max(double(CAL_RAD(:,:,3)),1)./(2^(16)-1)));

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
    for channel =1:3;
        button = questdlg(['Do you want to characterize the ' sprintf('%s',fct_channeltocolor(channel)) ' channel?'],'Mutichannel correction','Yes','No','Yes') ;
        if strcmp(button,'Yes')

            Iradtmp = Irad(:,:,channel);
            BCKGRNDtmp = BCKGRND(:,channel);
            
            button = questdlg('What kind of background analysis do you want to do?','Background analysis','Average','One-by-one','One-by-one') ;
            if strcmp(button,'Average')
                %ici ajouter un while au cas ou cela ne marche pas
                BCKGRNDtmp = BCKGRNDtmp*0+mean(BCKGRNDtmp);
            end
            
            DOSE = DOSE(:);
            [DOSE,isort] = sort(DOSE);
            RAWOD = zeros(nbfilms,1);
            
            for i=1:nbfilms
                [nlines,ncols] = size(Iradtmp);
                [xgrid,ygrid] = fct_gridindextopos(nlines,ncols,RES_RAD);
                xindex = fct_postoindex(xrad(i)-width_rect/2,xgrid);
                yindex = fct_postoindex(yrad(i)-height_rect/2,ygrid);
                A = imcrop(Iradtmp,[xindex yindex width_rect/RES_RAD height_rect/RES_RAD]);
                %A = fct_getroi(Iradtmp,RES_RAD,'Rectangular',[width_rect/RES_RAD height_rect/RES_RAD]);
                %figure; imagesc(A); set(gca,'DataAspectRatio',[1 1 1]);
                [mu, sigmamu, s] = fct_analyze_region(A);
                RAWOD(i) = mu;
            end
            
            OD = RAWOD(:)-BCKGRNDtmp(:);
            OD = OD(isort);
            
            %%%backup for autoconsitency test
            ODinit = OD;
            DOSEinit = DOSE;
            Npix = width_rect*height_rect/handles.CCDres/handles.CCDres;
            
            j=0;
            k = 0;
            option = 1:6;
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
            
%             %This code no longer works with 2014b (perhaps even before): #�"$%$&!!
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
            DOSEcal = DOSE;
            ODcal = OD;
            
            Nsave = N;
            optsave = opt;
            
            button='Yes';
            if strcmp(button,'Yes')
                showcurve = 0;
                n_vect = (handles.CCDres/min(height_rect,width_rect)).^2:5*(handles.CCDres/min(height_rect,width_rect)).^2:(handles.CCDres/0.1).^2;%it's a 1/N vector
                range = fliplr(handles.CCDres./sqrt(n_vect));
                sy_vect(1:length(range)) = 0;
                n_vect = (range/handles.CCDres).^2;
                hwait = waitbar(1/length(range),'Regression calculation');
                for j = 1:length(range)
                    RAWOD(1:nbfilms)=0;
                    for i=1:nbfilms
                        [nlines,ncols] = size(Iradtmp);
                        [xgridrad,ygridrad] = fct_gridindextopos(nlines,ncols,RES_RAD);
                        xindexrad = fct_postoindex(xrad(i)-range(j)/2,xgridrad);
                        yindexrad = fct_postoindex(yrad(i)-range(j)/2,ygridrad);
                        A_RAD = imcrop(Iradtmp,[xindexrad yindexrad range(j)/RES_RAD range(j)/RES_RAD]);
                        %A_RAD = fct_getroi(Iradtmp,RES_RAD,'Rectangular',[range(j)/RES_RAD range(j)/RES_RAD]);
                        [mu, sigmamu, s] = fct_analyze_region(A_RAD);
                        RAWOD(i) = mu;
                    end
                    OD = RAWOD-BCKGRNDtmp;
                    OD = OD(isort);
                    %%%%find min sy
                    SY = [];
                    option = 1:6;
                    Nrange = 2:min((nbfilms-2),10)
                    for opt = option
                        for N = Nrange
                            [p,sy,R,df] = fct_lsf(DOSE,OD,N,opt);
                            SY = [SY sy];
                        end
                    end
                    SY = SY';
                    %%%%
                    sy_vect(j) = min(SY);
                    if showcurve
                        h_tmp{j} = fct_show_calcurve(DOSE,OD,N,opt,0);
                    end
                    waitbar(j/length(range),hwait);
                end
                close(hwait);
                if showcurve
                    for j = 1:length(range)
                        if ishandle(h_tmp{j})
                            close(h_tmp{j});
                        end
                    end
                end
                x = 1./n_vect(:)+1./NBCK;
                y = sy_vect(:).^2;
                h_fit = figure('NumberTitle','off','Name','Uncertainty parameterization');
                p = polyfit(x,y,1);
                if p(1) < 0
                    p = [0 mean(y)];
                end
                sig0 = sqrt(p(2)/2);
                sig1 = sqrt(p(1));
                F = cat(2,x,x.^0);
                plot(x,F*p(:),'--b','Linewidth',2);
                hold on;
                errorbar(x,y,sy_vect.^2*sqrt(2/(length(DOSE)-N)),'ro','Linewidth',2);
                dmb = sy_vect.^2*sqrt(2/(length(DOSE)-N));
                [x(:) y(:) dmb(:)]
                dmb = F*p(:);
                [x(:) dmb(:)]
                hold off;
                dmb = sprintf('Linear fit y = 2(%.5f)^2 + (%.5f)^2x',sig0,sig1);
                legend({dmb,'Estimated variance'},'Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
                ylabel('Variance','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
                xlabel('1/N','Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
                set(gca,'YLim',[unique(min(0,min(sy_vect.^2)*1.35)) unique(max(sy_vect.^2)*1.35)]);
                set(gca,'XLim',[unique(min(0,min(1./n_vect+1./NBCK)*1.05)) unique(max(1./n_vect+1./NBCK)*1.05)]);
                set(gca,'Fontweight','Bold','Fontsize',12,'Fontname','Times New Roman');
                dmb = sprintf('L''incertitude est parametree sig0 = %.5f et sig1 = %.5f',sig0,sig1);
                h_msg = helpdlg(dmb,'Incertitude');
                uiwait(h_msg);
                close(h_fit);
                figure(handles.H);
            else
                sig0 = sycal/sqrt(2);
                sig1 = 0;
            end
            
            [ofilename,opathname]=uiputfile({'*.cal'},'Save calibration curve');
            size(DOSEcal)
            size(ODcal)
            
            if ~strcmp(class(ofilename),'double')
                file = fopen(fct_makecleanfilename(opathname,ofilename),'w');
                for i=1:nbfilms
                    fprintf(file,'%f\t%f\n',DOSEcal(i),ODcal(i));
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