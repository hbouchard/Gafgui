%function gamma = fct_GammaAnalysis(x,y,Z,x0,y0,Z0,distol,dostol);

%Must have resolution equal in both x and y for Z and Z0
%Z0 = reference image
%x0,y0 = positions vectors related to Z0
%Z = image to evaluate
%x,y = positions vectors related to Z
%distol,dostol: tolerances on distance and dose
%% Fake data for testing
clc;
clear;
close all;

res = 0.0169*10;
res0 = 0.05*10;

Xmax = 3;
Ymax = 5;
sig = 0.5;

distol = 0.3;
dostol = 0.03;

x = -Xmax:res:Xmax;
y = -Ymax:res:Ymax;
x0 = -Xmax:res0:Xmax;
y0 = -Ymax:res0:Ymax;

[X,Y] = meshgrid(x,y);
[X0,Y0] = meshgrid(x0,y0);

[lin,col] = size(X);
[lin0,col0] = size(X0);
xoffset = 0.5;%0.05*randn(lin,col);
yoffset = 0.5;%0.05*randn(lin,col);
Z0 = 1+exp(-0.5*((X0).^2+(Y0).^2)/sig^2) +...
     exp(-0.5*((X0).^2+(Y0-2).^2)/sig^2)+...
     exp(-0.5*((X0-1).^2+(Y0).^2)/sig^2);
maxval = max(Z0(:));
Z0 = Z0/maxval;
Z = 1+exp(-0.5*((X-xoffset).^2+(Y-yoffset).^2)/sig^2)+...
    exp(-0.5*((X-xoffset).^2+(Y-2-yoffset).^2)/sig^2)+...
    exp(-0.5*((X-1-xoffset).^2+(Y-yoffset).^2)/sig^2);
Z = Z/maxval;
Z = Z.*(1+0.01*randn(lin,col));

figure;
subplot(2,2,1);
imagesc(x,y,Z);
set(gca,'DataAspectRatio',[1 1 1]);
subplot(2,2,2);
imagesc(x0,y0,Z0);
set(gca,'DataAspectRatio',[1 1 1]);
% %this is fast as hell!
% subplot(1,2,2);
% xconv = conv2(X,ones(10,10));xconv = xconv(:,1);
% yconv = conv2(Y,ones(10,10));yconv = yconv(1,:);
% imagesc(xconv,yconv,conv2(Z,ones(10,10)));
% set(gca,'DataAspectRatio',[1 1 1]);
%% Basic info
%position information
x0min = min(x0); x0max = max(x0);
y0min = min(y0); y0max = max(y0);
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
lin0 = length(y0); col0 = length(x0);
lin = length(y); col = length(x);
resx = (xmax-xmin)/(col-1); resy = (ymax-ymin)/(lin-1);
resx0 = (x0max-x0min)/(col0-1); resy0 = (y0max-y0min)/(lin0-1);
%padding image with zeros: here assume the max x and y of both images is equal
if 0
    xpads = ceil((2*distol+resx)/resx);
    ypads = ceil((2*distol+resy)/resy);
    xmin = xmin-xpads*resx; xmax = xmax+xpads*resx;
    ymin = ymin-ypads*resy; ymax = ymax+ypads*resy;
    Z = cat(1,zeros(ypads,col),cat(1,Z,zeros(ypads,col)));
    Z = cat(2,zeros(lin+2*ypads,xpads),cat(2,Z,zeros(lin+2*ypads,xpads)));
    x = cat(2,xmin+(-1:-1:-xpads)*resx,cat(2,x,xmax+(1:xpads)*resx));
    y = cat(2,ymin+(-1:-1:-ypads)*resy,cat(2,y,ymax+(1:ypads)*resy));
    lin =length(y); col = length(x);
    subplot(2,2,1);
    imagesc(x,y,Z);
    set(gca,'DataAspectRatio',[1 1 1]);
    impixelinfo;
    [X,Y] = meshgrid(x,y);
    [lin,col] = size(X);
end
% convert position to matrix indeces
i = 1+round((y(:)-ymin)/resy); j = 1+round((x(:)'-xmin)/resx);
i0 = 1+round((y0(:)-y0min)/resy0); j0 = 1+round((x0(:)'-x0min)/resx0);
% convert matrix indeces to linear vector indeces
I = repmat(i(:),1,col); J = repmat(j(:)',lin,1);
K = I + (J-1)*lin;
I0 = repmat(i0(:),1,col0); J0 = repmat(j0(:)',lin0,1);
K0 = I0 + (J0-1)*lin0;
% reshape position/data matrices into vectors
X = reshape(X,lin*col,1);    Y = reshape(Y,lin*col,1);    Z = reshape(Z,lin*col,1);
X0 = reshape(X0,lin0*col0,1);Y0 = reshape(Y0,lin0*col0,1);Z0 = reshape(Z0,lin0*col0,1);
%

%%
tic
if 0
    Gamma = [];
    width = (2*ceil(2*distol/res))^2;
    %h = waitbar(0,'Please wait');
    for i=1:length(Z0)
        %waitbar(i/length(Z0),h);
        k = find((X0(i)-X).^2+(Y0(i)-Y).^2<=(2*distol)^2);
        Gammatmp = sqrt(((X0(i)-X(k)).^2+(Y0(i)-Y(k)).^2)/distol^2 + (Z0(i)-Z(k)).^2/dostol^2);
        tmp = 100*ones(1,width);
        tmp(1:length(Gammatmp)) = Gammatmp(:)';
        Gamma = cat(1,Gamma,tmp(:)');
    end
    %close(h);
elseif 0
    ZZ = repmat(Z(:)',length(Z0),1);
    XX = repmat(X(:)',length(X0),1);
    YY = repmat(Y(:)',length(Y0),1);
    ZZ0 = repmat(Z0(:),1,length(Z),1);
    XX0 = repmat(X0(:),1,length(X),1);
    YY0 = repmat(Y0(:),1,length(Y),1);  
    Gamma = sqrt(((XX-XX0).^2+(YY-YY0).^2)/distol^2+(ZZ-ZZ0).^2/dostol^2);
elseif 1
    for j=1:col
        ZZ = repmat(Z(:)',length(Z0),1);
        XX = repmat(X(:)',length(X0),1);
        YY = repmat(Y(:)',length(Y0),1);
        ZZ0 = repmat(Z0(:),1,length(Z),1);
        XX0 = repmat(X0(:),1,length(X),1);
        YY0 = repmat(Y0(:),1,length(Y),1);  
        Gamma = sqrt(((XX-XX0).^2+(YY-YY0).^2)/distol^2+(ZZ-ZZ0).^2/dostol^2);
    end
else
    Gamma = zeros(lin0*col0,1);
end
gamma = min(Gamma,[],2);
gamma = reshape(gamma,lin0,col0);
toc

subplot(2,2,3);
imagesc(x0,y0,gamma);
set(gca,'DataAspectRatio',[1 1 1]);
impixelinfo;
subplot(2,2,4);
imagesc(x0,y0,double(~logical(max(gamma-1,0))));
set(gca,'DataAspectRatio',[1 1 1]);
impixelinfo;