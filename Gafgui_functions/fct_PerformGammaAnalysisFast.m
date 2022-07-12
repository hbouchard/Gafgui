function gamma = fct_PerformGammaAnalysisFast(x,y,Z,x0,y0,Z0,distol,dostol)

[X,Y] = meshgrid(x,y);
[X0,Y0] = meshgrid(x0,y0);
[lin,col] = size(X);
[lin0,col0] = size(X0);

X = reshape(X,lin*col,1);    Y = reshape(Y,lin*col,1);    Z = reshape(Z,lin*col,1);
X0 = reshape(X0,lin0*col0,1);Y0 = reshape(Y0,lin0*col0,1);Z0 = reshape(Z0,lin0*col0,1);

ZZ = repmat(Z(:)',length(Z0(:)),1);
XX = repmat(X(:)',length(X0(:)),1);
YY = repmat(Y(:)',length(Y0(:)),1);
ZZ0 = repmat(Z0(:),1,length(Z(:)),1);
XX0 = repmat(X0(:),1,length(X(:)),1);
YY0 = repmat(Y0(:),1,length(Y(:)),1);  
Gamma = sqrt(((XX-XX0).^2+(YY-YY0).^2)/distol^2+(ZZ-ZZ0).^2/dostol^2);
gamma = min(Gamma,[],2);
gamma = reshape(gamma,lin0,col0);