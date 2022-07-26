% --------------------------------------------------------------------
function varargout = fct_GetCovarMatrixMulti(dose,Npixappl,filename,single,varargin)

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
ntappl = Fappl*p;

if Npixappl==0 %WB 2022: to be investigated why this could happen
    vappl = sycal^2; % voir comment choisir V1 et V2 en fonction de ce if la
else
    if single == 1
        if isempty(varargin) % absolute
            % this is for separate irradiations (rule 1)
            Vappl1 = diag(sycal^2+exp(polyval(qT,ntappl+t0)).^2./(Npixappl));
            % this is for same irradiation (rule 2)
            Vappl2 = sycal^2+diag(exp(polyval(qT,ntappl+t0)).^2./(Npixappl));
        elseif numel(varargin)==2 % relative
            norm = varargin{1};
            Npixnorm = varargin{2};

            tnorm = interp1(dose, ntappl, norm);
            Fnorm = fct_F_matrix(norm,M,type);
            dFnorm = fct_dF_matrix(norm,M,type);
            Fappl = cat(1,Fappl,Fnorm);
            dFappl = cat(1,dFappl, dFnorm);

            % This is for separate irradiations, as is the normalization (rule 1)
            Vappl1 = diag(sycal^2+exp(polyval(qT,ntappl+t0)).^2./(Npixappl));
            %extension
            Vappl1  = cat(2,cat(1,Vappl1 ,zeros(1,length(dose))),zeros(length(dose)+1,1));
            Vappl1(end,end)  = sycal^2  + exp(polyval(qT ,tnorm+t0)).^2./Npixnorm ;

            % This is for same irradiation, as is the normalization (rule 2)
            Vappl2 = sycal^2+diag(exp(polyval(qT,ntappl+t0)).^2./(Npixappl));
            %extension
            Vappl2  = cat(2,cat(1,Vappl2 ,sycal^2*ones(1,length(dose))) ,sycal^2*ones(length(dose)+1,1));
            Vappl2(end,end)  = sycal^2 + exp(polyval(qT ,tnorm+t0)).^2./(Npixnorm) ;
        end
        H = diag(1./(dFappl*p));
        varargout{1} = H* (Vappl1 + Fappl*Vp*Fappl') *H'; % rule 1
        varargout{2} = H* (Vappl2 + Fappl*Vp*Fappl') *H'; % rule 2
    elseif single == 2
        if isempty(varargin)  % absolute
           % Always separate irradiations so rule 1
                Vappl1 = diag(sycal^2+exp(polyval(qT,ntappl+t0)).^2/Npixappl);
        end
        H = diag(1./(dFappl*p));
        varargout{1} = H* (Vappl1 + Fappl*Vp*Fappl') *H'; % rule 1
    %ici je pense que j'ai implante comme il faut pour le cas ou on a
    %une dose de reference, alors il y a 2 valeurs de Npixappl
    %possibles, et dans ce cas 2 valeurs de dose. Ce serait bien de
    %faire un tableau avec les differents cas possibles ou on veut la
    %matrice de covariance, et pour chaque cas combien de valeurs de
    %Npixappl sont impliquees
    end
end
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