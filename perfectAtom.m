clear all;clc;
% This file generate a set of meta-atoms which covers entire 2 pi and has
% unity transmission.

fake = fopen('perfectAtom.txt','w');

r_list = [0.101:0.001:0.200];
Phase = linspace(0,2*pi,100);
index = linspace(1,100);
T=ones(1,100);
%for i = 1:100
%    T(i) = normpdf(index(i),100,20);
%end
for i=1:length(r_list)
    
    fprintf(fake,'%d ',r_list(i));
    fprintf(fake,'%d ',T(i));
    fprintf(fake,'%.4f\n',Phase(i));
    
end
fclose(fake);