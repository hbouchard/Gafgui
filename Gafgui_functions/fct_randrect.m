% --------------------------------------------------------------------
function x = fct_randrect(mu,a,halfopt)

if halfopt
    y = (0.5+rand(1,1)/2);
else
    y = rand(1,1);
end

% F = (x+a)/(2a);
x = 2*a*y-a;
x = x+mu;