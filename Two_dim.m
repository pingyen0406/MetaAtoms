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




clear all;
% Linux path
%folder_path = '/home/pingyen/Simulation/MATLAB/MetaAtoms/Lib562/60nmAl2O3/Al2O3_top/';
%addpath("/home/pingyen/Simulation/MATLAB/MetaAtoms/SubFunctions/");
% Windows path
folder_path = 'D:/Dropbox/MATLAB/MetaAtoms/Lib562/40nmAl2O3/';
addpath("D:/Dropbox/MATLAB/MetaAtoms/SubFunctions/");
fname_T = [folder_path,'SweepT562.txt'];
fname_Phase = [folder_path,'SweepPhase562.txt'];
outputlist = false;
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
f = 500; % focal length
beta =5; % beta angle of axicon(in degree)
size = 20; % radius or length of metalens (circle or square)
neff = 2.858; % effective index derived from FDTD
wavelength = 1.55;
decay = true;
decay_rate = 0.5;
x_range = [-100, 100];
x_res = 250;
y_range = [-100, 100];
y_res = 250;
z_range = [10, 1500];
z_res = 500;

%atomPos = squarePos("circle",[0,0],period,lens_radius);
atomPos = squarePos("square",[0,0],period,size);
N = sqrt(length(atomPos));
diff = [ones(1,N*N)*N*period;zeros(1,N*N)];
atomPos = cat(2,atomPos-diff,atomPos,atomPos+diff);
scatter(atomPos(1,:),atomPos(2,:));
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
focusing_field=Focusing_Slice(Dphase,T_list,atomPos,...
        x_range,z_range,0,x_res,z_res,1.55,true);
focal_field=Focal_Slice(Dphase,T_list,atomPos,...
    x_range,y_range,f,x_res,y_res,1.55,true); 
toc;
%}
%--------------------------------------------------------------------------



% Choosing desired part of phase data
Phase=NorPhase(Phase);
%%%%%%%%% Set an breakpoint this line to check the phase data.%%%%%%%

start_index=28;
stop_index=79;
Phase = Truncated_Phase(Phase,start_index,stop_index);
T = T(1,start_index:stop_index);
R_list = R_list(1,start_index:stop_index);

% Taking propagation phase into consideration 
%delay_phase = PropCorrect(length(atomPos),period,neff,wavelength);
%delay_phase = NorPhase(delay_phase);

tic;
% Creating focusing phase profile and doing interpolation
Dphase = SphericalOutput(0,f,[0,0],atomPos,1.55);
%Dphase = axiconOutput(0,beta,[0,0],atomPos,1.55);
[R_list,T_list]=Interpolation(Dphase,Phase,T,R_list);

% Simulating energy decay below the waveguide
if decay==true
    N_l = length(atomPos)/N;
    count=1;
    for i=1:N_l
        for j=1:N
            T_list(1,count)=(1-decay_rate*(i-1)/N_l)*T_list(1,count);
            count=count+1;
        end
        if count>N_l*N
            break
        end
    end
end


% Check the Radius and Transmission distribution
figure;
scatter3(atomPos(1,:),atomPos(2,:),T_list,'filled');
title("Transmission distribution");
figure;
scatter3(atomPos(1,:),atomPos(2,:),Dphase,'filled');
title("Phase distribution");


% Calculating the field and plot it out.(real, imag, and abs)
if plot_field==true
    focusing_field=Focusing_Slice(Dphase,T_list,atomPos,...
        x_range,z_range,0,x_res,z_res,1.55,false);
    focal_field=Focal_Slice(Dphase,T_list,atomPos,...
    x_range,y_range,f,x_res,y_res,1.55,false); 
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







