function gamma = fct_PerfromGammaAnalysisSlow(x,y,Z,res,x0,y0,Z0,res0,distol,dostol)

[X,Y] = meshgrid(x,y);
[X0,Y0] = meshgrid(x0,y0);
[lin,col] = size(X);
[lin0,col0] = size(X0);

X = reshape(X,lin*col,1);    Y = reshape(Y,lin*col,1);    Z = reshape(Z,lin*col,1);
X0 = reshape(X0,lin0*col0,1);Y0 = reshape(Y0,lin0*col0,1);Z0 = reshape(Z0,lin0*col0,1);

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
gamma = min(Gamma,[],2);
gamma = reshape(gamma,lin0,col0);