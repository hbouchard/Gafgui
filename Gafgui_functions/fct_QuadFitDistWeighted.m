function c = fct_QuadFitDistWeighted(x,y,i)

z = x(:)+y(:);
s = std(z);
mu = mean(z);
z = (z(:) - mu)/s;
p = 1/sqrt(2*pi)*exp(-0.5*z.^2);
w = sqrt(1./p(:));
w = min(w,sqrt(sqrt(2*pi)*exp(0.5*(3).^2)));
k = 1:length(z);
F = [w(k) w(k).*x(k) w(k).*x(k).^2];
q = F\(y(k).*w(k));
% diff = @(q) norm(y(k)-F*[abs(q(1)) q(2) q(3)]');
% options =  optimset('MaxFunEvals',10000,'MaxIter',10000,'Display','off');
% q = fminsearch(diff,[0 0 0],options);
c = q(i); 
