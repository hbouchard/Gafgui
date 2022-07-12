% --------------------------------------------------------------------
function [h,sycal,Errdose] = fct_ShowCalcurveMulti(DOSE,nTHETA,sTHETA,THETA0,Npix,res,channel,option,M,show_forcing)

h = -1;
syval = [];
Errdose = [];
if show_forcing==1
    Nrange =  2:min(max(size(DOSE))-1,10);
    optrange = 1:5;
else
    Nrange = M;
    optrange = option;
end

if show_forcing
    hwait = waitbar(1/length(optrange)/length(Nrange),'Fit calculation');
end
j = 0;
for opt = optrange
    for N = Nrange;
%         [p,sycal,R,df] = fct_lsf(DOSE,nTHETA,N,opt);
        [p,sycal,a1,a2,R,df] = fct_lsfmulti(DOSE,nTHETA,sTHETA,THETA0,N,opt);
        [validity,crit] = fct_iscurvevalid(p,N,opt,min(DOSE),max(DOSE));
        if validity||show_forcing
            %dosei = min(DOSE):1:max(DOSE);
            dosei = 1:(ceil(max(DOSE))-1);
            Fi = fct_F_matrix(dosei,N,opt);
            odi = Fi*p;
            dFi = fct_dF_matrix(dosei,N,opt);
            dodi = dFi*p;
            Q = Fi*inv(R);
            %maniere pratique de resoudre puisque la metode qr(F,0)
            %done une matrice R carree, on peut resoudre V = F*inv(R'*R)*F'
            %par V  = F*inv(R)*(F*inv(R))'=Q*Q'. Les elements diag sont
            %donc la lgine suivante
            sodi = sycal*sqrt(sum(Q.*Q,2));
            ErrnTHETA = mean(sodi);
            [Fcal,nametag] = fct_F_matrix(DOSE,N,opt);
            
            %%%%%%%%%
            %ici il ne s'agit pas du coefficient de
            %correlation!!!!!!!!!
            Rsq = sum((Fcal*p-mean(nTHETA)).^2)/sum((nTHETA-mean(nTHETA)).^2);
            %%%%%%%%%
            
            
            dFcal = fct_dF_matrix(DOSE,N,opt);
            
            K = length(DOSE);
            vect = dFcal*p ;
            H = zeros(K,K);
            for k = 1:K
                H(k,k) = 1/vect(k);
            end
            I = eye(K,K);
            V = H*(sycal^2*I+sycal^2*Fcal*inv(R'*R)*Fcal')*H;
            v = sqrt(diag(V));
            Errdose = mean(v);
            
            %%%%%%%
            if channel==1
                dmb = 'Red';
            elseif channel==2
                dmb = 'Green channel';
            elseif channel==3
                dmb = 'Blue channel';
            elseif channel==4
                dmb = 'Eigencolor ratio';
            end
            h = figure('NumberTitle','off','Name',[dmb ' calibration curve']);
            
            plot(DOSE,nTHETA,'ro','LineWidth',2);
            hold on;
            plot(dosei,odi,'--b','LineWidth',2);
            hold on;
            %plot(dosei,odi-sodi,'g--',DOSE(:)-v(:),Fcal*p,'r--',dosei,odi+sodi,'g--',DOSE(:)+v(:),Fcal*p,'r--','LineWidth',1);
            ylabel('Net optical density','Fontname','Times New Roman','Fontweight','bold','Fontsize',12);
            xlabel('Dose (CMU)','Fontname','Times New Roman','Fontweight','bold','Fontsize',12);
            set(gca,'YLim',[min(0,min(nTHETA)*1.05) max(nTHETA)*1.05]);
            set(gca,'XLim',[min(0,min(DOSE)*1.05) max(DOSE)*1.05]);
            str1 = '\theta';
            str2 = '\chi^2 per df';
            str3 = '\sigma_0';
            dumb=sprintf('Equation type: %s\nOrder: %d\nCalibration net %s uncertainty (%s): %s = %.4f\nAverage of estimated dose uncertainty = %.1f',nametag,N,str1,str2,str3,sycal,Errdose);
            %dumb=sprintf('Equation type: %s\nOrder: %d\nR^{2} = %.4f\nEstimated NnTHETA uncertainty (chi2) = %.4f\nMean estimated dose uncertainty = %.1f',nametag,N,Rsq,sycal,Errdose);
            %dumb=sprintf('Equation type = %s avec ordre %d\nR^{2} = %.4f\nSigma mesures estime (chi2) = %.3f',nametag,N,Rsq,sy);
            text(max(DOSE)*0.25,min(nTHETA)*0.1,dumb,'HorizontalAlignment','left','VerticalAlignment','bottom','Color','k','Fontname','Times New Roman','Fontweight','bold','Fontsize',12,'Edgecolor','k','Backgroundcolor','w');
            if ~validity
                dmb = sprintf('INVALIDE: 1 %d %d %d',crit(1),crit(2),crit(3));
                text((min(DOSE)+max(DOSE))/2,(min(nTHETA)+max(nTHETA))/2,dmb,'Color','k','Fontsize',26,'HorizontalAlignment','Center','FontWeight','Bold');
            end
            dumb1 = 'Uncertainty on estimated net \theta: 1 \sigma';
            dumb2 = 'Uncertainty on estimated 1 \sigma';
            %legend('Experimental data','Fitted data',dumb1,dumb2,4);
            legend('Experimental data','Fitted data','location','northwest');
            set(gca,'Fontname','Times New Roman','Fontweight','bold','Fontsize',12);
            hold off;
            grid on;
            [DOSE(:) nTHETA(:)]
            [dosei(:) odi(:)]
            
        end
        if show_forcing
            j = j+1;
            waitbar(j/length(optrange)/length(Nrange),hwait);
            figure(hwait);
        else
            %[DOSE(:) nTHETA(:)]
            %[dosei(:) odi(:)]
        end
    end
    
end
if show_forcing
    close(hwait);
end

h = figure('NumberTitle','off','Name',[dmb ' pixelwise standard deviation']);         
plot(nTHETA+THETA0,sTHETA,'ro','LineWidth',2);
hold on;
plot(nTHETA+THETA0,exp(a1+a2*(nTHETA+THETA0)),'--b','LineWidth',2);
hold on;
xlabel('Optical density','Fontname','Times New Roman','Fontweight','bold','Fontsize',12);
ylabel('Standard deviation (\sigma_1)','Fontname','Times New Roman','Fontweight','bold','Fontsize',12);
set(gca,'YLim',[min(sTHETA)*0.95 max(sTHETA)*1.05]);
set(gca,'XLim',[min(nTHETA+THETA0)*0.95 max(nTHETA+THETA0)*1.05]);

str = '\sigma^2 =\sigma_0^2+\sigma_1^2/N_{pix} ';

dumb = sprintf('Model: %s\nResolution = %f cm',str,res);
text(max(nTHETA+THETA0)*0.75,min(sTHETA),dumb,'HorizontalAlignment','left','VerticalAlignment','bottom','Color','k','Fontname','Times New Roman','Fontweight','bold','Fontsize',12,'Edgecolor','k','Backgroundcolor','w');

legend('Experimental data','Fitted data','location','northwest');
set(gca,'Fontname','Times New Roman','Fontweight','bold','Fontsize',12);
hold off;
grid on;