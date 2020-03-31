% --------------------------------------------------------------------
function [v,crit] = fct_iscurvevalid(p,N,type,dmin,dmax)

v = 1;
crit(1:3) = 0;
x = max(dmin,0.5):0.5:dmax;
%TEST1: function stricitly increasing
F = fct_F_matrix(x,N,type);
y = F*p;
test1 = norm(y - sort(y));
if test1==0
    v = v*1;
    crit(1) = 1;
else
    v = v*0;
end
%TEST2: 0 or 1 inflexion point
d2F = fct_d2F_matrix(x,N,type);
d2y = d2F*p;
s = sign(d2y);
test2 = s(1);
for n = 2:length(s)
    if s(n)~=s(n-1)
        test2 = cat(1,test2,s(n));
    end
end
%test2 = max(0,length(test2)-2); %allows 0 or 1 inflexion point
test2 = max(0,length(test2)-1); %allows 0 inflexion point
if test2==0
    v = v*1;
    crit(2) = 1;
else
    v = v*0;
end
%TEST3: no inflexion after y = 0.5*NODmax
k1 = find(abs(y-0.5)==min(abs(y-0.5)));
%k2 = find(abs(x-700)==min(abs(x-700)));
k2 = length(y);
s = d2y(k1:k2);
s = sign(s);
if(length(s)==0)
    v = v*0;%prevents that d2y=NaN
else
    test3 = s(1);
    for n = 2:length(s)
        if s(n)~=s(n-1)
            test3 = cat(1,test3,s(n));
        end
    end
    test3 = max(0,length(test3)-1);
    if test3==0
        v = v*1;
        crit(3) = 1;
    else
        v = v*0;
    end
end
if v==2
    F = fct_F_matrix(x,N,type);
    dF = fct_dF_matrix(x,N,type);
    d2F = fct_d2F_matrix(x,N,type);
    figure;
    plot(x,F*p);
    figure;
    plot(x,dF*p);
    figure;
    plot(x,d2F*p);
end