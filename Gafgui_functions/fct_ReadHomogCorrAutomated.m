% --------------------------------------------------------------------
 function [PR,PG,PB,x,rlims,glims,blims]  = fct_ReadHomogCorrAutomated(fname)

% clc;
% clear all;
% close all;
% 
% %
% fname = '/Users/Hugo/Desktop/junk.hm3';

file = fopen(fname,'r');
rlims = zeros(2,1);
glims = zeros(2,1);
blims = zeros(2,1);
rlims(1) = str2num(fgets(file));
rlims(2) = str2num(fgets(file));
glims(1) = str2num(fgets(file));
glims(2) = str2num(fgets(file));
blims(1) = str2num(fgets(file));
blims(2) = str2num(fgets(file));
nR = str2num(fgets(file));
nG = str2num(fgets(file));
nB = str2num(fgets(file));
PR = zeros(nR,1);
PG = zeros(nG,1);
PB = zeros(nB,1);
for n=1:nR
    PR(n) = str2num(fgets(file));
end
for n=1:nG
    PG(n) = str2num(fgets(file));
end
for n=1:nB
    PB(n) = str2num(fgets(file));
end
x = [];
line = 0;
while line~=-1
    line = fgets(file);
    if line~=-1
        x = cat(1,x,str2num(line));
    end
end
fclose(file);
rlims = sort(rlims);
glims = sort(glims);
blims = sort(blims);
