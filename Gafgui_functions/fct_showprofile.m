% --------------------------------------------------------------------
function f = fct_showprofile(x,profil,delta)

[i,j,k]=fct_caxindex(profil);
centerposition = x(j)
x = x-x(j);
%profil=profil/profil(j)*100;
s = fct_sym(profil);
h = fct_homog(profil);
FS=(k-i)*delta;
figure;
plot(x,profil,'LineWidth',2);

[x(:) profil(:)]

dmb = size(x);
dim = find(max(dmb)==dmb);
if dim==1
    xtmp=x;
else
    xtmp=x';
end

dmb=size(profil);
dim=find(max(dmb)==dmb);
if dim==1
    ytmp=profil;
else
    ytmp=profil';
end

A = cat(2,xtmp,ytmp);
fid=fopen('dummy.txt','w');
for ii=1:size(A,1)
    fprintf(fid,'%f %f\n',A(ii,1),A(ii,2));
end
fclose(fid);

hold on;
plot([x(i) x(k)],[profil(j)/2 profil(j)/2],'--k',[x(j) x(j)],[0 profil(j)],'--k','LineWidth',2);
hold on;
n=floor(0.8*(k-i+1)/2);
n=min(n,min(j-i,k-j));
plot([x(j-n) x(j-n)],[0 profil(j-n)],'--r',[x(j+n) x(j+n)],[0 profil(j+n)],'--r','LineWidth',1);
set(gca,'XLim',[min(x) max(x)]);

set(gca,'YLim',[0.9*min(profil) 1.1*max(profil)]);
xlabel('Off-axis distance (cm)');
ylabel('Value');
%grid on;
a=xlim;
b=ylim;
text(a(1),0.99*b(2),['Symmetry = ' num2str(s) '%'],'Color','k','Fontsize',12,'HorizontalAlignment','Left');
text(a(1),0.96*b(2),['Homogeneity = ' num2str(h) '%'],'Color','k','Fontsize',12,'HorizontalAlignment','Left');
text(a(1),0.93*b(2),['Field size = ' num2str(FS) ' cm'],'Color','k','Fontsize',12,'HorizontalAlignment','Left');
text(-0.2*(max(x)-min(x)),0.95*profil(j)/2,[num2str((x(j)-x(i))) ' cm'],'Color','k','Fontsize',12,'HorizontalAlignment','Center');
text(0.2*(max(x)-min(x)),0.95*profil(j)/2,[num2str((x(k)-x(j))) ' cm'],'Color','k','Fontsize',12,'HorizontalAlignment','Center');
f=1;
