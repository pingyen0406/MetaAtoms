clear all;clc;
% For Linux path
%folder_path = '/home/pingyen/Simulation/MATLAB/MetaAtoms/Lib562/60nmAl2O3/Al2O3_top/';
%addpath("/home/pingyen/Simulation/MATLAB/MetaAtoms/");
% For Windows path
folder_path = 'D:/Dropbox/MATLAB/MetaAtoms/Lib562/60nmAl2O3/Al2O3_top/';
addpath("D:/Dropbox/MATLAB/MetaAtoms/");
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
f = 46.8; % focal length
lens_radius = 50;
N = floor(2*lens_radius/period); % number of meta-atoms
neff = 2.858; % effective index derived from FDTD
wavelength = 1.55;
tmpPos = cell(N);
atomPos = zeros(2,1);
% crate a circular shape array of meta-atoms
for i=1:floor(N/2)
    for j=1:floor(N/2)
        if (period*i)^2+(period*j)^2 < lens_radius^2
            atomPos=cat(2,atomPos,[period*i;period*j]);
        else
            continue
        end               
    end
end 
atomPos=cat(2,atomPos,bsxfun(@times,atomPos,[-1;1]));
atomPos=cat(2,atomPos,bsxfun(@times,atomPos,[1;-1]));
    

% Choosing desired part of phase data
Phase=NorPhase(Phase);
% Set an breakpoint this line to check the phase data.
start_index=30;
stop_index=79;
Phase = Truncated_Phase(Phase,start_index,stop_index);
T = T(1,start_index:stop_index);
R_list = R_list(1,start_index:stop_index);

% Taking propagation phase into consideration 

%delay_phase = PropCorrect(N,lattice,neff,wavelength);
%delay_phase = NorPhase(delay_phase);

% Creating focusing phase profile and doing interpolation
Dphase = SphericalOutput(0,f,[0,0],atomPos,N,1.55);
[R_list,T_list]=Interpolation(delay_phase,Phase,T,R_list);

% Simulating energy decay below the waveguide
%{
alpha = 0.5;% decay rate 
for i=1:N
    T_list(1,i)=(1-alpha*i/N)*T_list(1,i);
end
%}

% Calculating the field and plot it out.(real, imag, and abs)
if plot_field==true
    Field=Focal_Slice(Dphase,T_list,atomPos,[-60,60],[-60,60],46.8,600,600,1.55);
end

% Output List
if outputlist==true 
    outf = fopen(outputname,'w');
    for i=1:length(R_list)
        fprintf(outf,'%f\n',R_list(i));
    end
    fclose(outf);
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


