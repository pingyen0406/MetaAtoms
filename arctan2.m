clear all;clc;clear workspace;
fname_T = 'SweepT800.txt';
fname_E = 'Enew.txt';
% (i,j) is 500+10*i nm height and 100+1*j nm radius
all_T = readmatrix(fname_T);
all_E = readmatrix(fname_E);
T = all_T;
E = all_E(:,1:end-1);
aPhase = zeros(size(E));
for i=1:size(E,1)
    aPhase(i,:) = atan2(imag(E(i,:)),real(E(i,:)));
    
end
H_range = [500:10:1490];
R_range = [30:2:248];
image(R_range,H_range,aPhase,'CDataMapping','scaled');
colormap('jet');

