clear all;
% For Linux path
%folder_path = '/home/pingyen/Simulation/MATLAB/MetaAtoms/Lib562/60nmAl2O3/Al2O3_top/';
%addpath("/home/pingyen/Simulation/MATLAB/MetaAtoms/SubFunctions/");
% For Windows path
folder_path = 'D:/Dropbox/MATLAB/MetaAtoms/Lib562/60nmAl2O3/Al2O3_top/';
addpath("D:/Dropbox/MATLAB/MetaAtoms/SubFunctions/");
fname_T = [folder_path,'SweepT562.txt'];
fname_Phase = [folder_path,'SweepPhase562.txt'];
outputlist = false;
outputname = [folder_path,'spherical1200.txt'];
plot_field = true;

% Read data from S4
% (i,j) is 500+10*i nm height and 100+1*j nm radius
all_T = readmatrix(fname_T);
all_Phase = readmatrix(fname_Phase);
height = 1200;
T = all_T(height/10-50+1,:);
Phase = all_Phase(height/10-50+1,:);
R_list = [0.03:0.002:0.248];

% Parameters
period = 0.562;
f = 100; % focal length
L = 40;
N = floor(L/period); % number of meta-atoms
neff = 2.858; % effective index derived from FDTD
wavelength = 1.55;
atomPos = zeros(2,N);
% crate x-position array of meta-atoms
for i=1:N
    x_now = i-floor(N/2)-1;
    atomPos(1,i)=period*x_now;
    
end 
%---------------------------- Test data-------------------------------
%{
test_name = 'perfectAtom.txt';
test_f = readmatrix(test_name);
test_r = transpose(test_f(:,1));
test_T = transpose(test_f(:,2));
test_Phase = transpose(test_f(:,3));
[sphericalList,Dphase] = SphericalOutput(test_Phase,test_T,test_r,f,lattice,N,1.55);
amp_list = interp1(test_r,test_T,sphericalList);
Field=Eatom(Dphase,amp_list,atomPos,[0,60],[1,20],0,1200,400,1.55);
%}
%-----------------------------------------------------------------------

tic;
% Choosing desired part of phase data
Phase=NorPhase(Phase);
% Set an breakpoint this line to check the phase data.
start_index=1;
stop_index=77;
Phase = Truncated_Phase(Phase,start_index,stop_index);
T = T(1,start_index:stop_index);
R_list = R_list(1,start_index:stop_index);

% Taking propagation phase into consideration 
%delay_phase = PropCorrect(N,lattice,neff,wavelength);
%delay_phase = NorPhase(delay_phase);

% Creating focusing phase profile and doing interpolation
Dphase = SphericalOutput(0,f,[0,0],atomPos,1.55);
[R_list,T_list]=Interpolation(Dphase,Phase,T,R_list);

% Simulating energy decay below the waveguide
%{
alpha = 0.5;% decay rate 
for i=1:N
    T_list(1,i)=(1-alpha*i/N)*T_list(1,i);
end
%}

% Calculating the field and plot it out.(real, imag, and abs)
if plot_field==true
    focal_field=Focal_Slice(Dphase,T_list,atomPos,[-30,30],[-30,30],f,300,300,1.55);
    focusing_field=Focusing_Slice(Dphase,T_list,atomPos,[-30,30],[1,151],0,300,750,1.55);
end

% Output List
if outputlist==true 
    outf = fopen(outputname,'w');
    for i=1:length(R_list)
        fprintf(outf,'%f\n',R_list(i));
    end
    fclose(outf);
end
toc;
% Check the Radius and Transmission distribution
figure;
scatter3(atomPos(1,:),atomPos(2,:),T_list,'filled');
title("Transmission distribution");


