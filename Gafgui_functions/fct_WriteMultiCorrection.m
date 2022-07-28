% --------------------------------------------------------------------
function err = fct_WriteMultiCorrection(fname,Rot,rrange,grange,brange,trange)
% function err = fct_WriteMultiCorrection(fname,Rot,rrange,grange,brange,trange,znorm)
err = 1;
%
file = fopen(fname,'w');
fprintf(file,'%f %f %f\n',Rot(1,1),Rot(1,2),Rot(1,3));
fprintf(file,'%f %f %f\n',Rot(2,1),Rot(2,2),Rot(2,3));
fprintf(file,'%f %f %f\n',Rot(3,1),Rot(3,2),Rot(3,3));
fprintf(file,'%d %d %d\n',min(rrange),max(rrange),0);
fprintf(file,'%d %d %d\n',min(grange),max(grange),0);
fprintf(file,'%d %d %d\n',min(brange),max(brange),0);
fprintf(file,'%f %f %d',min(trange),max(trange),0);
% fprintf(file,'%f %f %d',min(trange),max(trange),znorm);
fclose(file);
%
err = 0;