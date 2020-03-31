% --------------------------------------------------------------------
function err = fct_WriteCorrmatrixAuto(fname,xx,s1,s2,s3,c1,c2,c3,p1,p2,p3,order1,order2,sym,m,n,bits)

err = 1;
file = fopen(fname,'w');
fprintf(file,'%d\n',bits);
fprintf(file,'%d\n',order1);
fprintf(file,'%d\n',order2);
fprintf(file,'%d\n',sym);
fprintf(file,'%d\n',m);
fprintf(file,'%d\n',n);
fprintf(file,'%d\n',length(p1));
for i=1:length(p1)
    fprintf(file,'%e\n',p1(i));
end
for i=1:length(p2)
    fprintf(file,'%e\n',p2(i));
end
for i=1:length(p3)
    fprintf(file,'%e\n',p3(i));
end
fprintf(file,'%d\n',length(xx));
for i=1:length(xx)
    fprintf(file,'%e\n',xx(i));
end
for i=1:length(s1)
    fprintf(file,'%e\n',s1(i));
end
for i=1:length(s2)
    fprintf(file,'%e\n',s2(i));
end
for i=1:length(s3)
    fprintf(file,'%e\n',s3(i));
end
for i=1:length(c1)
    fprintf(file,'%e\n',c1(i));
end
for i=1:length(c2)
    fprintf(file,'%e\n',c2(i));
end
for i=1:length(c3)-1
    fprintf(file,'%e\n',c3(i));
end
fprintf(file,'%e',c3(end));
fclose(file);
err = 0;