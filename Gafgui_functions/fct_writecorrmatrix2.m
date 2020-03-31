% --------------------------------------------------------------------
function err = fct_writecorrmatrix2(fname,x,Y)
err = 1;
N = length(Y);
file = fopen(fname,'w');
fprintf(file,'2\n');
fprintf(file,'%d\n',N);
fprintf(file,'%f\n',x);
for i=1:N
    y = Y{i};
    fprintf(file,'%f\n',y);
end
fclose(file);
err = 0;