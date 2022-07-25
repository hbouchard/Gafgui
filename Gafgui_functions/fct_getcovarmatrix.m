% --------------------------------------------------------------------
function [V1, V2] = fct_getcovarmatrix(dose,Npixappl,Npixbck,filename,single,varargin)

%here Npixappl can be of length 1 or same length of dose
%length Npixbck must be 1

file = fopen(filename,'r');
% WB 2022: need here to use fct_READCAlFileMulti and not fct_readcalfile
%[DOSE,OD,M,type,sigparam,dummy,t0] = fct_readcalfile(file);
[DOSE,OD,s1,t0,dummy,res,channel,type,M] = fct_ReadCalFileMulti(file);
fclose(file);
[p,sycal,R,df] = fct_lsf(DOSE,OD,M,type);
qT  = polyfit(OD+t0, log(s1) ,1);

Fappl = fct_F_matrix(dose,M,type);
dFappl = fct_dF_matrix(dose,M,type);
Fcal = fct_F_matrix(DOSE,M,type);
%N = length(find(DOSE~=0));
N = length(DOSE);
Vp = (inv(Fcal'*Fcal)*Fcal')*(sycal^2*eye(N,N))*(inv(Fcal'*Fcal)*Fcal')';
tappl = Fappl*p;

if Npixappl==0
    vappl = sycal^2; % voir comment choisir V1 et V2 en fonction de ce if la
else
    l = length(Npixappl);
    if l ==1
        if single == 1
            if isempty(varargin) % absolute
                Vappl1 = 2*sycal^2+diag(exp(polyval(qT,tappl+t0)).^2./(Npixappl+Npixbck));
                Vappl2 = diag(2*sycal^2+exp(polyval(qT,tappl+t0)).^2./(Npixappl+Npixbck));
            elseif length(varargin)==2 % relative
                norm = varargin{1};
                Npixnorm = varargin{2};

                tnorm = interp1(dose, tappl, norm);
                Fnorm = fct_F_matrix(norm,M,type);
                dFnorm = fct_dF_matrix(norm,M,type);
                Fappl = cat(1,Fappl,Fnorm);
                dFappl = cat(1,dFappl, dFnorm);

                Vappl1 = 2*sycal^2+diag(exp(polyval(qT,tappl+t0)).^2./(Npixappl+Npixbck));
                %extension
                Vappl1  = cat(2,cat(1,Vappl1 ,2*sycal^2*ones(1,length(dose))) ,2*sycal^2*ones(length(dose)+1,1));
                Vappl1(end,end)  = 2*sycal^2 + 2*exp(polyval(qT ,tnorm+t0)).^2./(Npixnorm) ;

                Vappl2 = diag(2*sycal^2+exp(polyval(qT,tappl+t0)).^2./(Npixappl+Npixbck));
                %extension
                Vappl2  = cat(2,cat(1,Vappl2 ,zeros(1,length(dose))),zeros(length(dose)+1,1));
                Vappl2(end,end)  = 2*sycal^2  + 2*exp(polyval(qT ,tnorm+t0)).^2./Npixnorm ;
                
            end
        %elseif
        end
    elseif l==2 % a compléter si jamais
        vappl = 2*sycal^2+s1.^2.*(1/Npixappl(1)+1/Npixbck);
        vappl = [vappl 2*sycal^2+s1.^2.*(1/Npixappl(2)+1/Npixbck)];
    end
    %ici je pense que j'ai implante comme il faut pour le cas ou on a
    %une dose de reference, alors il y a 2 valeurs de Npixappl
    %possibles, et dans ce cas 2 valeurs de dose. Ce serait bien de
    %faire un tableau avec les differents cas possibles ou on veut la
    %matrice de covariance, et pour chaque cas combien de valeurs de
    %Npixappl sont impliquees
end

H = diag(1./(dFappl*p));
V1 = H* (Vappl1 + Fappl*Vp*Fappl') *H'; % rule 1
V2 = H* (Vappl2 + Fappl*Vp*Fappl') *H'; % rule 2

end


% if length(vappl)~=length(dose)
%     vappl = vappl*ones(length(dose),1);
% end
% n = length(dose);
% Vappl = zeros(n,n);
% for i=1:n
%     Vappl(i,i) = vappl(i);
% end
%
% Fappl = fct_F_matrix(dose,M,type);
% dFappl = fct_dF_matrix(dose,M,type);
% Fcal = fct_F_matrix(DOSE,M,type);
% N = length(dose);
% vect = dFappl*p ;
% H = zeros(N,N);
% for k = 1:N
%     H(k,k) = 1/vect(k);
% end
% I = eye(N,N);
% V = H*(Vappl+sycal^2*Fappl*inv(R'*R)*Fappl')*H;