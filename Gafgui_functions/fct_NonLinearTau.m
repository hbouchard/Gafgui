function [tau,m,n] = fct_NonLinearTau(xmat,a,b,c)

tau = [];
for j=1:numel(xmat)
    x = -log(a+b*xmat{j} +c*xmat{j}.^2);
    tau = cat(1,tau,x/mean(mean(x)));
end
[m,n] = size(tau);
tau = tau(:);