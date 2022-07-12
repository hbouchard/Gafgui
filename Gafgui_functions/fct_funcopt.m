function val = fct_funcopt(theta,t,D,X)

val = abs(theta - X(1).*log(t+X(4)) - X(2).*log(t+X(4)).*((D+X(5)).*log(D+X(5))-(D+X(5))) - X(4).*D(:) );