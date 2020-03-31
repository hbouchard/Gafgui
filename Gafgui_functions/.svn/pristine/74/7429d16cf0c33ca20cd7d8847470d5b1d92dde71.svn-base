

% --------------------------------------------------------------------
function dF = fct_dF_matrix(x,M,type)

dF = [];
A = 100;
B = 500;
alpha = 0.01;

if type==1
    
    for n = 1:M
        df = n*x.^(n-1);
        dF = cat(2,dF,df(:));
    end
    
elseif type==2
    df = 1./(x+A);
    dF = cat(2,dF,df(:));
    for n = 2:M
        df = n*(log(x/A+1).^(n-1))./(x+A);
        dF = cat(2,dF,df(:));
    end
    
elseif type==3
    f = B./(x.^2+B.^2);
    dF = cat(2,dF,f(:));
    for n = 2:M
        f = n.*atan(x/B).^(n-1).*B./(x.^2+B.^2);
        dF = cat(2,dF,f(:));
    end
    
elseif type==4
    df = x.*0+1;
    dF = cat(2,dF,df(:));
    for n = 1:M
        df = ((n+3)/2)*x.^((n+1)/2);
        dF = cat(2,dF,df(:));
    end
    
elseif type==5
    
    for m = 1:M
        S1 = [];
        S2 = [];
        for n = 0:m-1
            S1 = cat(2,S1,(alpha*x(:)).^n/factorial(n));
        end
        if m>1
            for n = 0:m-2
                S2 = cat(2,S2,(alpha*x(:)).^n/factorial(n));
            end
        else
            S2 = x(:)*0;
        end
        S1 = sum(S1,2);
        S2 = sum(S2,2);
        df = alpha*exp(-alpha*x(:)).*S1 - alpha*exp(-alpha*x(:)).*S2 ;
        dF = cat(2,dF,df(:));
    end
    
elseif type==6
    
    df = x.*0+1;
    dF = cat(2,dF,df(:));
    df = M*x.^(M-1);
    dF = cat(2,dF,df(:));
    
end

dF = fliplr(dF);