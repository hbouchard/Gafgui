

% --------------------------------------------------------------------
function d2F = fct_d2F_matrix(x,M,type)

d2F = [];
A = 100;
B = 500;
alpha = 0.01;

if type==1
    d2f = 0*x;
    d2F = cat(2,d2F,d2f(:));
    for n = 2:M
        d2f = n*(n-1)*x.^(n-2);
        d2F = cat(2,d2F,d2f(:));
    end
    
elseif type==2
    d2f = -1./((x+A).^2);
    d2F = cat(2,d2F,d2f(:));
    for n = 2:M
        d2f = (n-1)*n*(log(x/A+1).^(n-2))./((x+A).^2);
        d2f = d2f - n*(log(x/A+1).^(n-1))./((x+A).^2);
        d2F = cat(2,d2F,d2f(:));
    end
    
elseif type==3
    f = -B./((x.^2+B.^2).^2).*2.*x;
    d2F = cat(2,d2F,f(:));
    for n = 2:M
        f = B.^2*(n-1)*n./(x.^2+B.^2).^2.*atan(x/B).^(n-2) - B*2*x.*n./(x.^2+B.^2).^2.*atan(x/B).^(n-1);
        d2F = cat(2,d2F,f(:));
    end
    
elseif type==4
    d2f = x.*0;
    d2F = cat(2,d2F,d2f(:));
    for n = 1:M
        d2f = ((n+3)/2)*((n+1)/2)*x.^((n-1)/2);
        d2F = cat(2,d2F,d2f(:));
    end
    
elseif type==5
    
    for m = 1:M
        S1 = [];
        S2 = [];
        S3 = [];
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
        if m>2
            for n = 0:m-3
                S3 = cat(2,S3,(alpha*x(:)).^n/factorial(n));
            end
        else
            S3 = x(:)*0;
        end
        S1 = sum(S1,2);
        S2 = sum(S2,2);
        S3 = sum(S3,2);
        d2f = -alpha^2*exp(-alpha*x(:)).*S1 + 2*alpha^2*exp(-alpha*x(:)).*S2 - alpha^2*exp(-alpha*x(:)).*S3 ;
        d2F = cat(2,d2F,d2f(:));
    end
    
elseif type==6
    
    d2f = x.*0;
    d2F = cat(2,d2F,d2f(:));
    if M>1
        d2f = (M-1)*M*x.^(M-2);
    else
        d2f = x(:)*0;
    end
    d2F = cat(2,d2F,d2f(:));
    
end

d2F = fliplr(d2F);