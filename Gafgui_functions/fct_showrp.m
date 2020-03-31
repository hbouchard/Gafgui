% --------------------------------------------------------------------
function f = fct_showrp(x,rp,delta)

rp=rp/max(rp)*100;
isurf= find(x==0);
imax = find(rp==max(rp));
RP=rp;
if imax>1
    rp(1:imax-1)=rp(imax);
end
ininety=min(find(abs(rp-90)==min(abs(rp-90))));
ieigtyfive=min(find(abs(rp-85)==min(abs(rp-85))));
ieigty=min(find(abs(rp-80)==min(abs(rp-80))));
ififty=min(find(abs(rp-50)==min(abs(rp-50))));
isurf= min(find((abs(x)==min(abs(x)))));
imax = max(find(rp==max(rp)));
iten=min(find(abs(x-10)==min(abs(x-10))));
itwenty=min(find(abs(x-20)==min(abs(x-20))));
Rmax=delta*(imax-isurf);

figure;
plot(x,RP,'LineWidth',2);
set(gca,'XLim',[min(x) max(x)]);
set(gca,'YLim',[0 1.1*max(RP)]);
xlabel('Off-axis distance (cm)');
ylabel('Dose (%)');
a=xlim;
b=ylim;
xlabel('Depth (cm)');
ylabel('Dose (%)');
hold on;

Dsurf=RP(isurf);
text(x(isurf),RP(isurf),sprintf('Dsurf = %.1f %%',Dsurf),'Color','k','Fontsize',12,'HorizontalAlignment','Left');
plot([x(isurf) x(isurf)],[b(1) RP(isurf)],'--k',[a(1) x(isurf)],[RP(isurf) RP(isurf)],'--k','LineWidth',1);
hold on;

text(x(imax),RP(imax),sprintf('R100 = %.2f cm',Rmax),'Color','k','Fontsize',12,'HorizontalAlignment','Left');
plot([x(imax) x(imax)],[b(1) RP(imax)],'--k',[a(1) x(imax)],[RP(imax) RP(imax)],'--k','LineWidth',1);
hold on;

if min(rp)<=20
    
    if min(rp)<=90
        Rninety=delta*(ininety-isurf);
        text(x(ininety),RP(ininety),sprintf('R90 = %.2f cm',Rninety),'Color','k','Fontsize',12,'HorizontalAlignment','Left');
        plot([x(ininety) x(ininety)],[b(1) RP(ininety)],'--k',[a(1) x(ininety)],[RP(ininety) RP(ininety)],'--k','LineWidth',1);
        hold on;
    end
    if min(rp)<=85
        Reigtyfive=delta*(ieigtyfive-isurf);
        text(x(ieigtyfive),RP(ieigtyfive),sprintf('R85 = %.2f cm',Reigtyfive),'Color','k','Fontsize',12,'HorizontalAlignment','Left');
        plot([x(ieigtyfive) x(ieigtyfive)],[b(1) RP(ieigtyfive)],'--k',[a(1) x(ieigtyfive)],[RP(ieigtyfive) RP(ieigtyfive)],'--k','LineWidth',1);
        hold on;
    end
    if min(rp)<=80
        Reigty=delta*(ieigty-isurf);
        text(x(ieigty),RP(ieigty),sprintf('R80 = %.2f cm',Reigty),'Color','k','Fontsize',12,'HorizontalAlignment','Left');
        plot([x(ieigty) x(ieigty)],[b(1) RP(ieigty)],'--k',[a(1) x(ieigty)],[RP(ieigty) RP(ieigty)],'--k','LineWidth',1);
        hold on;
    end
    if min(rp)<=50
        Rfifty=delta*(ififty-isurf)
        ififty
        size(rp)
        min(rp)
        text(x(ififty),RP(ififty),sprintf('R50 = %.2f cm',Rfifty),'Color','k','Fontsize',12,'HorizontalAlignment','Left');
        plot([x(ififty) x(ififty)],[b(1) RP(ififty)],'--k',[a(1) x(ififty)],[RP(ififty) RP(ififty)],'--k','LineWidth',1);
        hold on;
    end
else
    if max(x)>=10
        Dten=RP(iten);
        text(x(iten),RP(iten),sprintf('D10 = %.1f %%',Dten),'Color','k','Fontsize',12,'HorizontalAlignment','Left');
        plot([x(iten) x(iten)],[b(1) RP(iten)],'--k',[a(1) x(iten)],[RP(iten) RP(iten)],'--k','LineWidth',1);
        hold on;
    end
    if max(x)>=20
        Dtwenty=RP(itwenty);
        text(x(itwenty),RP(itwenty),sprintf('D10 = %.1f %%',Dtwenty),'Color','k','Fontsize',12,'HorizontalAlignment','Left');
        plot([x(itwenty) x(itwenty)],[b(1) RP(itwenty)],'--k',[a(1) x(itwenty)],[RP(itwenty) RP(itwenty)],'--k','LineWidth',1);
        hold on;
    end
end
