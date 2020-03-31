function [c,sigmax] = fct_LinearFitDistWeighted(x,y,i)

z = x(:)+y(:);
s = std(z);
mu = mean(z);
z = (z(:) - mu)/s;
% p = 1/sqrt(2*pi)*exp(-0.5*z.^2);
% w = sqrt(1./p(:));
w = sqrt(sqrt(2*pi)*exp(0.5*z.^2));
N = length(z);
%erfgauss(x) = erf(x/sqrt(2)) = P(|z|<x);
%If Gaussian, probability shouldn't go beyond 1-1/N
sigmax = sqrt(2)*erfinv(1-1/N);
%We won't allow sigmax to be over a certain value to avoid badly conditioned system
%I am not sure what is the critical value, but pmin ~ 1e-10 yields 5 orders
%of magnitudes differences in weights
pmin = 1/sqrt(2*pi)*1e-10;%1e-15 is near the limit so 1-pmin ~= 1;
sigcritical = sqrt(2)*erfinv(1-pmin);
sigmax = min(sigmax,sigcritical);
lim = sqrt(sqrt(2*pi)*exp(0.5*sigmax.^2));
w = min(w,lim);
%We treat all values above lim as aberrant
w = (1-double(logical(z>sigmax))).*w;
k = find(w~=0);
F = [w(k) w(k).*x(k)];
q = F\(y(k).*w(k));

s = std(y(k));
yhat = [x(k).^0 x(k).^1]*q;
chi2perdf = sum(w(k).^2.*(yhat-y(k)).^2/s^2)/sum(w(k).^2);
q = cat(1,q,chi2perdf);
%for i=1, it is the y intersect, for i=2, it is the slope
c = q(i); 
