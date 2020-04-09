% Calculating the scalar field emitted by 2D meta-atom array. 
% General using method:
% 1.Reading the phase and transmission data from simulation(S4).
% 2.Setting the parameters. (period, lens size, focal length...)
% 3.Creating 2xN matrix that contain position of every meta-atom.
% 4.Using NorPhase.m to check the phase data. And using Truncated_Phase.m
%   to pick out the interval if needed.
% 5.Using SphericalOutput.m or axiconOutput.m to generate the desired
%   phase profile of every meta-atom. And then using Interploation.m to 
%   find the corresponding radius and transmission.
% 6.Using Focal_slice.m or Focusing_slice.m to calaulate the field at given
%   position. If there is symmetric property of the lens, the computation
%   time can reduce to 1/2 or 1/4 time.




clear all; close all;
% Linux path
%folder_path = '/home/pingyen/Simulation/MATLAB/MetaAtoms/Lib562/60nmAl2O3/Al2O3_top/';
%addpath("/home/pingyen/Simulation/MATLAB/MetaAtoms/SubFunctions/");
% Windows path
folder_path = 'D:/Dropbox/MATLAB/MetaAtoms/Lib562/60nmAl2O3/Al2O3_top/';
addpath("D:/Dropbox/MATLAB/MetaAtoms/SubFunctions/");
fname_T = [folder_path,'SweepT562.txt'];
fname_Phase = [folder_path,'SweepPhase562.txt'];
outputlist = true;
outputname = [folder_path,'axicon1200.txt'];
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
f = 10000; % focal length
beta =31; % beta angle of axicon(in degree)
lens_radius = 28.1; % radius or length of metalens
N = floor(2*lens_radius/period); % number of meta-atoms
neff = 2.858; % effective index derived from FDTD
wavelength = 1.55;
atomPos = zeros(0,0);
tmpPos = zeros(0,0);

atomPos = squarePos("circle",period,lens_radius);
atomPos2 = squarePos("square",period,3);

figure;
scatter(atomPos(1,:),atomPos(2,:),'filled');
figure;
scatter(atomPos2(1,:),atomPos2(2,:),'filled');
%---------------------------- Test data-----------------------------------
%{
tic;
test_name = 'perfectAtom.txt';
test_f = readmatrix(test_name);
test_r = transpose(test_f(:,1));
test_T = transpose(test_f(:,2));
test_Phase = NorPhase(transpose(test_f(:,3)));
test_Phase(1,end)=0;
Dphase = SphericalOutput(0,f,[0,0],atomPos,1.55);
[R_list,T_list]=Interpolation(Dphase,test_Phase,test_T,test_r);
focal_field=Focal_Slice(Dphase,T_list,atomPos,[-25,25],[-25,25],f,100,100,1.55);
focusing_field=Focusing_Slice(Dphase,T_list,atomPos,[-25,25],[1,201],0,100,2000,1.55);
toc;
%}
%--------------------------------------------------------------------------



% Choosing desired part of phase data
Phase=NorPhase(Phase);
%%%%%%%%% Set an breakpoint this line to check the phase data.%%%%%%%

start_index=35;
stop_index=79;
Phase = Truncated_Phase(Phase,start_index,stop_index);
T = T(1,start_index:stop_index);
R_list = R_list(1,start_index:stop_index);

% Taking propagation phase into consideration 
delay_phase = PropCorrect(length(atomPos),period,neff,wavelength);
delay_phase = NorPhase(delay_phase);
tic;
% Creating focusing phase profile and doing interpolation
Dphase = axiconOutput(delay_phase,beta,[0,0],atomPos,1.55);
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
    focusing_field=Focusing_Slice(Dphase,T_list,atomPos,...
        [-30,30],[1,101],0,600,1000,1.55,true);
    %focal_field=Focal_Slice(Dphase,T_list,atomPos,...
    %[-200,200],[-200,200],f,1000,1000,1.55,true); 
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
%figure;
%scatter3(atomPos(1,:),atomPos(2,:),T_list,'filled');
%title("Transmission distribution");
%figure;
%scatter3(atomPos(1,:),atomPos(2,:),R_list,'filled');
%title("Transmission distribution");



