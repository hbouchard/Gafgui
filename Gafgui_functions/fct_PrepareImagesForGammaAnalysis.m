function [x,y,Z] = fct_PrepareImagesForGammaAnalysis(x0,y0,x,y,Z,distol,option);

%position information
x0min = min(x0); x0max = max(x0);
y0min = min(y0); y0max = max(y0);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
lin0 = length(y0); col0 = length(x0);
lin = length(y); col = length(x);
resx = (xmax-xmin)/(col-1); resy = (ymax-ymin)/(lin-1);
resx0 = (x0max-x0min)/(col0-1); resy0 = (y0max-y0min)/(lin0-1);
%padding image with zeros
if option
    delx = max(0,max(x0min-2*distol-xmin,x0max+2*distol-xmax));
    dely = max(0,max(y0min-2*distol-ymin,y0max+2*distol-ymax));
    xpads = ceil((delx+resx)/resx);
    ypads = ceil((dely+resy)/resy);
    xmin = xmin-xpads*resx; xmax = xmax+xpads*resx;
    ymin = ymin-ypads*resy; ymax = ymax+ypads*resy;
    Z = cat(1,zeros(ypads,col),cat(1,Z,zeros(ypads,col)));
    Z = cat(2,zeros(lin+2*ypads,xpads),cat(2,Z,zeros(lin+2*ypads,xpads)));
    x = cat(2,xmin+(-1:-1:-xpads)*resx,cat(2,x,xmax+(1:xpads)*resx));
    y = cat(2,ymin+(-1:-1:-ypads)*resy,cat(2,y,ymax+(1:ypads)*resy));
    lin =length(y); col = length(x);
end


