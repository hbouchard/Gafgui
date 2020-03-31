% --------------------------------------------------------------------
function err = fct_WriteCorrmatrixSemiAuto(fname,xi,ri,gi,bi,r0i,g0i,b0i,X,R,G,B,R0,G0,B0)

err = 1;
file = fopen(fname,'w');
[m,n] = size(xi);
fprintf(file,'%.d\n',m);
%1
    for j=1:n; fprintf(file,'%.1f ',xi(1,j));end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',xi(i,j)); end;
end
%2
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',ri(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',ri(i,j)); end;
end
%3
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',gi(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',gi(i,j)); end;
end
%4
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',bi(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',bi(i,j)); end;
end
%5
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',r0i(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',r0i(i,j)); end;
end
%6
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',g0i(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',g0i(i,j)); end;
end
%7
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',b0i(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',b0i(i,j)); end;
end
%reset
fprintf(file,'\n');
[m,n] = size(X);
fprintf(file,'%.d\n',m);
%1
    for j=1:n; fprintf(file,'%.1f ',X(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',X(i,j)); end;
end
%2
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',R(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',R(i,j)); end;
end
%3
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',G(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',G(i,j)); end;
end
%4
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',B(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',B(i,j)); end;
end
%5
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',R0(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',R0(i,j)); end;
end
%6
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',G0(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',G0(i,j)); end;
end
%7
fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',B0(1,j)); end;
for i=2:m
    fprintf(file,'\n');
    for j=1:n; fprintf(file,'%.1f ',B0(i,j)); end;    
end
err = 0;