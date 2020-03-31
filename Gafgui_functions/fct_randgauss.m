% --------------------------------------------------------------------
function x = fct_randgauss(mu,sigma,halfopt)

if halfopt
    y = (0.5+rand(1,1)/2);
else
    y = rand(1,1);
end

y0 = 0;
del = 10;
x0 = 0;
s = 1;
while abs(y0-y)>1e-6
    del = del/10;
    x = x0-del*10:del:x0+del*10;
    F = erf(x/(sqrt(2)))/2+0.5;
    k = find(abs(F-y)==min(abs(F-y)));
    x0 = x(k);
    y0 = erf(x0/(sqrt(2)*s))/2+0.5;
end

x = x0*sigma + mu;