
% --------------------------------------------------------------------
function [F,nametag] = fct_F_matrix(x,M,type)

F = [];
A = 100;
B = 500;
alpha = 0.01;

if type==1
    
    for n = 1:M
        f = x.^n;
        F = cat(2,F,f(:));
    end
    nametag = 'y = \Sigma a_nx^n';
    
elseif type==2
    
    for n = 1:M
        f = log(x/A+1).^n;
        F = cat(2,F,f(:));
    end
    nametag = 'y = \Sigma a_nln(x/A+1)^n';
    
elseif type==3
    
    for n = 1:M
        f = atan(x/B).^n;
        F = cat(2,F,f(:));
    end
    nametag = 'y = \Sigma a_natan(x/B)^n';
    
elseif type==4
    f = x;
    F = cat(2,F,f(:));
    for n = 1:M
        f = x.^((n+3)/2);
        F = cat(2,F,f(:));
    end
    nametag = 'y = a_0x + \Sigma a_nx^(^n^+^3^)^/^2';
    
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
        f = 1-exp(-alpha*x(:)).*S1;
        F = cat(2,F,f(:));
    end
    nametag = 'Weighted multiple hit theory';
    
elseif type==6
    
    f = x;
    F = cat(2,F,f(:));
    f = x.^M;
    F = cat(2,F,f(:));
    nametag = 'y = a_1x + a_2x^n';
    
end

F = fliplr(F);