function err = fct_WriteMultichCorrLinearNew(fname,ISOD,DOSERANGE,TAURANGE,pr,pg,pb,b,res)
err = 1;
file = fopen(fname,'w')
fprintf(file,'%e\n',ISOD);
fprintf(file,'%e\n',min(DOSERANGE));
fprintf(file,'%e\n',max(DOSERANGE));
fprintf(file,'%e\n',min(TAURANGE));
fprintf(file,'%e\n',max(TAURANGE));
fprintf(file,'%e\n',res);
fprintf(file,'%e\n',b);
fprintf(file,'%e\n',length(pr));
for i = 1:length(pr)
    fprintf(file,'%e\n',pr(i));
end
for i = 1:length(pg)
    fprintf(file,'%e\n',pg(i));
end
for i = 1:length(pb)
    fprintf(file,'%e\n',pb(i));
end
fclose(file);
err = 0;