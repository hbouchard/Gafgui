function [S,p] = fct_NonLinearRatioStd(X,Y,Z,q,a,b)

% p = (1+atan(q)*2/pi)/20;
    p = q;
    p(1) = min(max(q(1),-abs(a)),abs(a));
    p(3) = min(max(q(3),-abs(a)),abs(a));
    p(5) = min(max(q(5),-abs(a)),abs(a));
    p(2) = min(max(q(2),-abs(b)),0.0);
    p(4) = min(max(q(4),-abs(b)),0.0);
    p(6) = min(max(q(6),-abs(b)),0.0);
    

r1 = (X-log10(1+p(1)*10.^(-X))+p(2))./(Z-log10(1+p(5)*10.^(-Z))+p(6));
r2 = (Y-log10(1+p(3)*10.^(-Y))+p(4))./(Z-log10(1+p(5)*10.^(-Z))+p(6));
S = sum(std(r1,[],1).^2./mean(r1,1).^2) + sum(std(r2,[],1).^2./mean(r2,1).^2);
% 
% ax = (1+atan(q(1))*2/pi)/2;
% bx = (1+atan(q(2))*2/pi)/2;
% ay = (1+atan(q(3))*2/pi)/2;
% by = (1+atan(q(4))*2/pi)/2;
% az = (1+atan(q(5))*2/pi)/2;
% bz = (1+atan(q(6))*2/pi)/2;
%
% s = zeros(length(DOSE),1);
% for j=1:length(DOSE)
%     k=find(D==DOSE(j));
%     r1 = (X-log10(1+ax*exp(-X))+bx)./(Z-log10(1+az*exp(-Z))+bz);
%     r2 = (Y-log10(1+ay*exp(-Y))+by)./(Z-log10(1+az*exp(-Z))+bz);
%     s(j) = std(r1)/mean(r1)+std(r1)/mean(r1);
% end
% S = sum(s);
% 
