% --------------------------------------------------------------------
function V = fct_getcovarmatrix(dose,Npixappl,Npixbck,filename)

%here Npixappl can be of length 1 or same length of dose
%length Npixbck must be 1

file = fopen(filename,'r');
[DOSE,OD,M,type,sigparam] = fct_readcalfile(file);
fclose(file);
[p,sycal,R,df] = fct_lsf(DOSE,OD,M,type);
if Npixappl==0
    vappl = sycal^2;
else
    l = length(Npixappl);
    if l ==1
        vappl = 2*sigparam(1)^2+sigparam(2)^2.*(1/Npixappl+1/Npixbck);
    elseif l==2
        vappl = 2*sigparam(1)^2+sigparam(2)^2.*(1/Npixappl(1)+1/Npixbck);
        vappl = [vappl 2*sigparam(1)^2+sigparam(2)^2.*(1/Npixappl(2)+1/Npixbck)];
    end
    %ici je pense que j'ai implante comme il faut pour le cas ou on a
    %une dose de reference, alors il y a 2 valeurs de Npixappl
    %possibles, et dans ce cas 2 valeurs de dose. Ce serait bien de
    %faire un tableau avec les differents cas possibles ou on veut la
    %matrice de covariance, et pour chaque cas combien de valeurs de
    %Npixappl sont impliquees
end
if length(vappl)~=length(dose)
    vappl = vappl*ones(length(dose),1);
end
n = length(dose);
Vappl = zeros(n,n);
for i=1:n
    Vappl(i,i) = vappl(i);
end
Fappl = fct_F_matrix(dose,M,type);
dFappl = fct_dF_matrix(dose,M,type);
Fcal = fct_F_matrix(DOSE,M,type);
N = length(dose);
vect = dFappl*p ;
H = zeros(N,N);
for k = 1:N
    H(k,k) = 1/vect(k);
end
I = eye(N,N);
V = H*(Vappl+sycal^2*Fappl*inv(R'*R)*Fappl')*H;