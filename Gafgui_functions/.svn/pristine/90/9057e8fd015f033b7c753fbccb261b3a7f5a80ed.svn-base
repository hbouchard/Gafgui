
% --------------------------------------------------------------------
function [optimum,inflex] = fct_polyanalysis(p,domain)

order = length(p)-1;
n = order:-1:0;
dp = p(:).*n(:);
dp = dp(1:order);
n = (order-1):-1:0;
d2p = dp(:).*n(:);
d2p = d2p(1:order-1);
r = roots(dp); %roots of derivatve
k = find(imag(r)==0);
optimum = [];
for i=1:length(k)
    if (r(k(i))<=max(domain))&&(r(k(i))>=min(domain))
        optimum = [optimum r(k(i))];
    end
end

r = roots(d2p); %roots of 2nd derivatve
k = find(imag(r)==0);
inflex = [];
for i=1:length(k)
    if (r(k(i))<=max(domain))&&(r(k(i))>=min(domain))
        inflex = [inflex r(k(i))];
    end
end
%5 AOUT
%HERE IMPLEMENT: if the function finishes in values at edge that are
%greater than central maximum, it is wrong