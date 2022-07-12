function gamma = fct_GammaAnalysis2D(x,y,Z,x0,y0,Z0,distol,dostol,maxdatasizeMB,maxsearchdis);

%Z0 = reference image
%x0,y0 = positions vectors related to Z0
%Z = image to evaluate
%x,y = positions vectors related to Z
%distol,dostol: tolerances on distance and dose

%% Calculation
datasize = length(Z0(:))*length(Z(:))*16;%in Bytes
%maxdatasizeMB = 64e6;
[lin0,col0] = size(Z0);
gamma = zeros(lin0,col0);
h = waitbar(0,'Please wait');
%
ndiv = ceil(sqrt(datasize/(maxdatasizeMB*1e6)));
imin = 0; imax = 0;  
for i=1:ndiv
    waitbar(i/ndiv,h);
    imin = min(lin0,imax+1); 
    imax = min(lin0,imin+ceil(lin0/ndiv)-1);
    %restricted ROI is a square not circle
    ii = intersect(find(y<=(y0(imax)+maxsearchdis)),find(y>=(y0(imin)-maxsearchdis)));
    jmin = 0; jmax = 0;
    for j=1:ndiv
        jmin = min(col0,jmax+1); 
        jmax = min(col0,jmin+ceil(col0/ndiv)-1);
        if (imax-imin)*(jmax-jmin)
            %restricted ROI is a square not circle
            jj = intersect(find(x<=(x0(jmax)+maxsearchdis)),find(x>=(x0(jmin)-maxsearchdis)));
            x0red = x0(jmin:jmax);
            y0red = y0(imin:imax);
            Z0red = Z0(imin:imax,jmin:jmax);
            xred = x(jj);
            yred = y(ii);
            Zred = Z(ii,jj);
            gammared = fct_PerformGammaAnalysisFast(xred,yred,Zred,x0red,y0red,Z0red,distol,dostol);
            gamma(imin:imax,jmin:jmax) = gammared;
        end
    end
end
close(h);