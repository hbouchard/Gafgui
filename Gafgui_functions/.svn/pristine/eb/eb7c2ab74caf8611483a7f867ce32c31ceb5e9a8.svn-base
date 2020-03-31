% --------------------------------------------------------------------
function h = fct_homog(y)

[i,j,k]=fct_caxindex(y);
n=floor(0.8*(k-i+1)/2);
n=min(n,min(j-i,k-j));

if size(y,1)>size(y,2)
    y=y';
end

hm(1:2*n+1)=y(j-n:j+n);
h=((max(hm)-min(hm))/2/y(j)+1)*100;