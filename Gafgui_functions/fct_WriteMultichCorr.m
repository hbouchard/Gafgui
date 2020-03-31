function err = fct_WriteMultichCorr(fname,ISOD,DOSERANGE,TAURANGE,pr,pg,pb,pwr,pwg,pwb)
err = 1;
file = fopen(fname,'w')
fprintf(file,'%e\n',ISOD);
fprintf(file,'%e\n',min(DOSERANGE));
fprintf(file,'%e\n',max(DOSERANGE));
fprintf(file,'%e\n',min(TAURANGE));
fprintf(file,'%e\n',max(TAURANGE));
fprintf(file,'%e\n',length(pr));
fprintf(file,'%e\n',length(pwr));
for i = 1:length(pr)
    fprintf(file,'%e\n',pr(i));
end
for i = 1:length(pg)
    fprintf(file,'%e\n',pg(i));
end
for i = 1:length(pb)
    fprintf(file,'%e\n',pb(i));
end
for i = 1:length(pwr)
    fprintf(file,'%e\n',pwr(i));
end
for i = 1:length(pwg)
    fprintf(file,'%e\n',pwg(i));
end
for i = 1:length(pwb)
    fprintf(file,'%e\n',pwb(i));
end
fclose(file);
err = 0;