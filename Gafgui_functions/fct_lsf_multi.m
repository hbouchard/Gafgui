% --------------------------------------------------------------------
function [a,Va] = fct_lsf_multi(DOSE,OD,sOD,N,type,isnoabscissa)

x = DOSE(:);
y = OD(:);
w = 1./sOD;

%Create function matrix and do LSF
F = fct_F_matrix2(x,N,type,isnoabscissa);
[n,m] = size(F);
W = [];
for i=1:m
    W = cat(2,W,w(:));
end
M = inv((F.*W)'*(F.*W))*(F.*W)';
a = M*(y.*w);
Va =  M*M';
r = y - F*a;
