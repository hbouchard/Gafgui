% --------------------------------------------------------------------
function str = fct_meanval_output(DOSE,V)
mu = DOSE(1);
smu = sqrt(V(1,1));
str = sprintf('Absolute value = %.1f CMU\n----------------',mu);
str = sprintf('%s\nUncertainties:\ntype A = +/- %.1f (%.1f %%)',str,smu,100*smu/mu);
%Ici les source d'incertitude sont causee par le setup et l'output de
%la machine. On a environ 1.7% pour l'output (dist. rect. +/-3%) et on
%a l'erreur de setup +/- 1 mm qui varie avec le gradient utilis�. On
%consid�re environ 1.5% mais cela peut-�tre plus
sB = sqrt(1.7^2 + 1.5^2)/100;
% str = sprintf('%s\ntype B (approx.) = +/- %.1f (%.1f %%)',str,sB*mu,100*sB);
% str = sprintf('%s\ntotal (approx.) = +/- %.1f (%.1f %%)',str,sqrt(smu^2 + (sB*mu)^2),100*sqrt((smu/mu)^2 + sB^2));
if max(size(DOSE))~=1
    r = DOSE(1)/DOSE(2);
    sr = r*sqrt(V(1,1)/DOSE(1)^2+V(2,2)/DOSE(2)^2 - 2*V(1,2)/DOSE(1)/DOSE(2));
    str = sprintf('%s\n\n\nRelative value = %.4f\n----------------',str,r);
    str=sprintf('%s\nUncertainties (k=1):\ntype A = +/- %.4f (%.1f %%)',str,sr,sr/r*100);
%     str = sprintf('%s\ntype B (approx.) = +/- %.4f (%.1f %%)',str,sB*r,100*sB);
%     str = sprintf('%s\ntotal (approx.) = +/- %.4f (%.1f %%)',str,sqrt(sr^2 + (sB*r)^2),100*sqrt((sr/r)^2 + sB^2));
    
end