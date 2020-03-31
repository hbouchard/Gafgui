
% --------------------------------------------------------------------
function x = fct_randtriang(mu,a,halfopt)

if halfopt
    y = (0.5+rand(1,1)/2);
else
    y = rand(1,1);
end

if y==0.5
    x = 0;
elseif y<0.5
    x = sqrt(2*a*a*y)-a; %F = (x+a)^2/(2*a^2) ---> x = sqrt(2a^2F)-a
else
    x = a*(1-sqrt(2-2*y)); %F = -(x^2-2ax-a^2)/(2*a^2) ---> x = a(1-sqrt(2(1-F)))
end

x = x+mu;