
% --------------------------------------------------------------------
function err = fct_writecorrmatrix3(fname,xlim,x0,slim,p,npol,orders)

err = 1;
file = fopen(fname,'w');
fprintf(file,'3\n');
fprintf(file,'%e\n',min(xlim));
fprintf(file,'%e\n',max(xlim));
fprintf(file,'%e\n',x0);
fprintf(file,'%e\n',min(slim));
fprintf(file,'%e\n',max(slim));
fprintf(file,'%d',npol);
for j=1:npol
    fprintf(file,'\n%d',orders(j));
end
for j=1:npol
    tmp = p{j};
    for i = 1:length(tmp)
        fprintf(file,'\n%e',tmp(i));
    end
end
fclose(file);
err = 0;