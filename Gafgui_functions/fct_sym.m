% --------------------------------------------------------------------
function s = fct_sym(y)

[i,j,k]=fct_caxindex(y);
n=floor(0.8*(k-i+1)/2);
n=min(n,min(j-i,k-j));

if size(y,1)>size(y,2)
    y=y';
end

u(1:2*n+1)=y(j-n:j+n);
a=u./fliplr(u);
b=fliplr(u)./u;
sm=max(a(:),b(:));
s=max(sm)*100;