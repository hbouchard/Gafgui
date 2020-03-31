% --------------------------------------------------------------------
function x = fct_RandStandGaussian(N)

P = rand(N,1);
%P(x) = 0.5+0.5*erf(x/sqrt(2))
x = erfinv(2*P-1)*sqrt(2);
