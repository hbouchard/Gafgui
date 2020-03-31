% --------------------------------------------------------------------
function err = fct_writecorrmatrix(fname,orders,xmin,xmax,smin,smax,p)

err = 1;
npol = length(orders);
file = fopen(fname,'w');
fprintf(file,'1\n');
fprintf(file,'%e\n',npol);
for i = 1:npol
    fprintf(file,'%e\n',orders(i));
end
fprintf(file,'%e\n',xmin);
fprintf(file,'%e\n',xmax);
fprintf(file,'%e\n',smin);
fprintf(file,'%e\n',smax);
for i = 1:npol
    tmp = p{i};
    tmp = tmp(:)';
    fprintf(file,'%e\n',tmp);
end
fclose(file);
err = 0;