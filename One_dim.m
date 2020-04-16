clear all;
% For Linux path
%folder_path = '/home/pingyen/Simulation/MATLAB/MetaAtoms/Lib562/60nmAl2O3/Al2O3_top/';
%addpath("/home/pingyen/Simulation/MATLAB/MetaAtoms/SubFunctions/");
% For Windows path
folder_path = 'D:/Dropbox/MATLAB/MetaAtoms/Lib562/40nmAl2O3/';
addpath("D:/Dropbox/MATLAB/MetaAtoms/SubFunctions/");
fname_T = [folder_path,'SweepT562.txt'];
fname_Phase = [folder_path,'SweepPhase562.txt'];
outputlist = false;
outputname = [folder_path,'spherical1200_prop.txt'];
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
beta = 5; % beta angle for axicon lens
f = 2000; % focal length
L = 56.2; % array length (um)
N = floor(L/period)+1; % number of meta-atoms
neff = 2.858; % effective index derived from FDTD
wavelength = 1.55;
atomPos = zeros(2,N);
decay =true;
decay_rate = 0.5;
x_range = [-75, 75];
x_res = 250;
y_range = [-75, 75];
y_res = 250;
z_range = [10, 1500];
z_res = 500;



%---------------------------- Test data-------------------------------
%{
test_name = 'perfectAtom.txt';
test_f = readmatrix(test_name);
test_r = transpose(test_f(:,1));
test_T = transpose(test_f(:,2));
test_Phase = transpose(test_f(:,3));
tt = NorPhase(test_Phase);
[sphericalList,Dphase] = SphericalOutput(test_Phase,test_T,test_r,f,lattice,N,1.55);
amp_list = interp1(test_r,test_T,sphericalList);
Field=Eatom(Dphase,amp_list,atomPos,[0,60],[1,20],0,1200,400,1.55);
%}
%-----------------------------------------------------------------------

tic;
% Choosing desired part of phase data
Phase=NorPhase(Phase);
% Set an breakpoint this line to check the phase data.
start_index=35;
stop_index=79;
Phase = Truncated_Phase(Phase,start_index,stop_index);
T = T(1,start_index:stop_index);
R_list = R_list(1,start_index:stop_index);

% crate 1D position array of meta-atoms
for i=1:N
    if mod(N,2)~=0
        x_now = i-floor(N/2)-1;
        atomPos(1,i)=period*x_now;
    else
        x_now = i-floor(N/2);
        atomPos(1,i)=period*x_now;
        atomPos(1,i) = atomPos(1,i)-0.5*period;
    end
    
end

% Taking propagation phase into consideration
%delay_phase = PropCorrect(length(atomPos),period,neff,wavelength);
%delay_phase = NorPhase(delay_phase);

tic;



% Creating focusing phase profile and doing interpolation
%Dphase = SphericalOutput(0,f,[0,0],atomPos,1.55);
Dphase = axiconOutput(0,beta,[0,0],atomPos,1.55);
[R_list,T_list]=Interpolation(Dphase,Phase,T,R_list);

% Simulating energy decay below the waveguide
if decay==true 
    for i=1:N
        T_list(1,i)=(1-decay_rate*(i-1)/N)*T_list(1,i);
    end
end
%}

% Check the Transmission and phase distribution
figure;
scatter(atomPos(1,:),T_list,'filled');
title("Transmission distribution");
figure;
scatter(atomPos(1,:),Dphase,'filled');
title("Phase distribution");

% Calculating the field and plot it out.(real, imag, and abs)
if plot_field==true
    focusing_field=Focusing_Slice(Dphase,T_list,atomPos,...
        x_range,z_range,0,x_res,z_res,1.55,false);
    %focal_field=Focal_Slice(Dphase,T_list,atomPos,...
    %x_range,y_range,f,x_res,y_res,1.55,false); 
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



