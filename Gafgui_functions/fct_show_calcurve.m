% --------------------------------------------------------------------
function [h,sycal,Errdose] = fct_show_calcurve(DOSE,OD,M,option,show_forcing)

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
        [p,sycal,R,df] = fct_lsf(DOSE,OD,N,opt);
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
            ErrOD = mean(sodi);
            [Fcal,nametag] = fct_F_matrix(DOSE,N,opt);
            
            %%%%%%%%%
            %ici il ne s'agit pas du coefficient de
            %correlation!!!!!!!!!
            Rsq = sum((Fcal*p-mean(OD)).^2)/sum((OD-mean(OD)).^2);
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
            
            h = figure('NumberTitle','off','Name','Calibration curve');
            
            plot(DOSE,OD,'ro','LineWidth',2);
            hold on;
            plot(dosei,odi,'--b','LineWidth',2);
            hold on;
            %plot(dosei,odi-sodi,'g--',DOSE(:)-v(:),Fcal*p,'r--',dosei,odi+sodi,'g--',DOSE(:)+v(:),Fcal*p,'r--','LineWidth',1);
            ylabel('Net optical density','Fontname','Times New Roman','Fontweight','bold','Fontsize',12);
            xlabel('Dose (CMU)','Fontname','Times New Roman','Fontweight','bold','Fontsize',12);
            set(gca,'YLim',[min(0,min(OD)*1.05) max(OD)*1.05]);
            set(gca,'XLim',[min(0,min(DOSE)*1.05) max(DOSE)*1.05]);
            dumb=sprintf('Equation type: %s\nOrder: %d\nCalibration OD uncertainty (chi2) = %.4f\nMean estimated dose uncertainty = %.1f',nametag,N,sycal,Errdose);
            %dumb=sprintf('Equation type: %s\nOrder: %d\nR^{2} = %.4f\nEstimated NOD uncertainty (chi2) = %.4f\nMean estimated dose uncertainty = %.1f',nametag,N,Rsq,sycal,Errdose);
            %dumb=sprintf('Equation type = %s avec ordre %d\nR^{2} = %.4f\nSigma mesures estime (chi2) = %.3f',nametag,N,Rsq,sy);
            text(max(DOSE)*0.25,min(OD)*0.1,dumb,'HorizontalAlignment','left','VerticalAlignment','bottom','Color','k','Fontname','Times New Roman','Fontweight','bold','Fontsize',12,'Edgecolor','k','Backgroundcolor','w');
            if ~validity
                dmb = sprintf('INVALIDE: 1 %d %d %d',crit(1),crit(2),crit(3));
                text((min(DOSE)+max(DOSE))/2,(min(OD)+max(OD))/2,dmb,'Color','k','Fontsize',26,'HorizontalAlignment','Center','FontWeight','Bold');
            end
            dumb1 = 'Uncertainty on estimated NOD: 1 \sigma';
            dumb2 = 'Uncertainty on estimated 1 \sigma';
            %legend('Experimental data','Fitted data',dumb1,dumb2,4);
            legend('Experimental data','Fitted data','location','northwest');
            set(gca,'Fontname','Times New Roman','Fontweight','bold','Fontsize',12);
            hold off;
            grid on;
            [DOSE(:) OD(:)]
            [dosei(:) odi(:)]
            
        end
        if show_forcing
            j = j+1;
            waitbar(j/length(optrange)/length(Nrange),hwait);
            figure(hwait);
        else
            %[DOSE(:) OD(:)]
            %[dosei(:) odi(:)]
        end
    end
    
end
if show_forcing
    close(hwait);
end