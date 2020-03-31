% function handles = fct_CreateSignalToODMap(handles);

clc;
clear all;
close all;

fct_AddGafguiFctPath();
handles = fct_initGafgui();

button = questdlg('Do you want to start from the beginning?','Shortcut','Yes','No','Yes') ;
% figure(handles.H);

if strcmp(button,'No')
    flag = 1;
else
    flag = 2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%Get data
if flag >= 2
    
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
        button = questdlg('How do you wish to enter OD values?','OD values','File','Manual','File') ;
        flag = 1;
        OD = [];
        if strcmp(button,'Manual')
            for j = 1:nprompt
                clear prompt;
                for i= 1+5*(j-1):min(5*j,nbfilms)
                    dumb = sprintf('Film #%d\nOD',i);
                    prompt{i-5*(j-1)} = dumb;
                end
                answer = inputdlg(prompt,'OD values',1);
                OD = cat(1,OD,str2double(answer));
            end
            OD = OD';
            %ici il devrait y avoir une facon de valider les valeurs de dose
            OD = abs(OD);
        else
            [ifilename,ipathname] = uigetfile({'*.txt'},'Choose file containing OD values');
            if ~strcmp(class(ifilename),'double')
                file = fopen(fct_makecleanfilename(ipathname,ifilename),'r');
                OD = fscanf(file,'%f',[1 inf]);
            else
                flag = 0;
            end
        end
        if flag~=0
            button = questdlg('Do you want to save this step?','Shortcut','Yes','No','Yes') ;
            if strcmp(button,'Yes')
                %ici ajouter un while au cas ou cela ne marche pas
                [ofilename,opathname]=uiputfile({'*.mat'},'Work to save');
                if ofilename==0
                else
                    filename = fct_makecleanfilename(opathname,ofilename);
                    save(filename,'posx','posy','rmat','bmat','gmat','bits','OD');
                end
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
if flag ~=-1
    nbbands = numel(rmat);
    r = []; g = [];  b = [];
    for j=1:nbbands
       r = cat(1,r,mean(rmat{j}(:))/(2^bits-1));
       g = cat(1,g,mean(gmat{j}(:))/(2^bits-1));
       b = cat(1,b,mean(bmat{j}(:))/(2^bits-1));
    end
    [r,k] = sort(r);  g = g(k);  b = b(k);
    OD = sort(OD(:));
    I = sort(10.^(-OD(1:nbbands)));    
    r1 = polyval(polyfit(I(1:2),r(1:2),1),0); r2 = polyval(polyfit(I((end-1):end),r((end-1):end),1),1,'linear');
    g1 = polyval(polyfit(I(1:2),g(1:2),1),0); g2 = polyval(polyfit(I((end-1):end),g((end-1):end),1),1,'linear');
    b1 = polyval(polyfit(I(1:2),b(1:2),1),0); b2 = polyval(polyfit(I((end-1):end),b((end-1):end),1),1,'linear');
    r = [r1 r(:)' r2]';
    g = [g1 g(:)' g2]';
    b = [b1 b(:)' b2]';
    I = [0  I(:)' 1]';
    %%%
    figure;
    plot(I,r,'.r',I,g,'.g',I,b,'.b');
    legend('Red','Green','Blue','location','northwest');
    xlabel('Relative light');
    ylabel('Relative scanner signal');
    set(gca,'DataAspectRatio',[1 1 1]);
    %%%
    [ofilename,opathname]=uiputfile({'*.od'},'Save OD calibration curve');
    
    if ~strcmp(class(ofilename),'double')
        file = fopen(fct_makecleanfilename(opathname,ofilename),'w');
        for i=1:(nbbands-1)
            fprintf(file,'%f\t%f\t%f\t%f\n',I(i),r(i),g(i),b(i));
        end
        fprintf(file,'%f\t%f\t%f\t%f\n',I(nbbands),r(nbbands),g(nbbands),b(nbbands));
        fclose(file);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%