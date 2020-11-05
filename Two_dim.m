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
input=jsondecode(fileread('input.json'));
% Linux path
%folder_path = '/home/pingyen/Simulation/MATLAB/MetaAtoms/Lib562/60nmAl2O3/Al2O3_top/';
%addpath("/home/pingyen/Simulation/MATLAB/MetaAtoms/SubFunctions/");
% Windows path
folder_path = input.folder_path;
current_fol = pwd;
addpath([current_fol,'\','SubFunctions\']);
fname_T = [folder_path,input.inputName_T];
fname_Phase = [folder_path,input.inputName_Phase];
plot_focalField = input.plot_focalField;
plot_focusingField = input.plot_focusingField;
outputlist = input.outputlist;
outputname = [input.output_path,input.output_name];

% Read data from S4 calculation
% (i,j) is 500+10*i nm height and 100+1*j nm radius
all_T = readmatrix(fname_T);
all_Phase = readmatrix(fname_Phase);
height = input.height;
T = all_T(height/0.01-50+1,:);
Phase = all_Phase(height/0.01-50+1,:);
R_range = [0.05:0.002:0.348];

% Parameters
period = input.period;
lens_type = input.lens_spec.lens_type; % "axicon" or "spherical"
lens_size = input.size;
% 'Radius' or 'Length' of metalens (circle or square)
f_num = 3.5; % f-number
f = input.lens_spec.f; % focal length
gt_angle = input.lens_spec.gt_angle; % grating angle
center=input.center; % middle point of the pattern
beta =input.lens_spec.beta; % beta angle of axicon(in degree)
neff = input.neff; % effective index derived from FDTD
wavelength = 1.55;
phase_delay=input.phase_delay;
decay = input.dacay;
decay_rate = 0.5; % The ramaining power in the waveguide. 

% phase shift between [0,1], you can modify this to have a better intensity
% distribution.
phi0 = 0.4; 

[atomPos_X,atomPos_Y] = squarePos([0,0],period,lens_size);
size_array = size(atomPos_X);

%diff = [ones(1,N*N)*N*period;zeros(1,N*N)];
% Use 3 squares to create a rectangle(20um*60um lens)
%atomPos = cat(2,atomPos-diff,atomPos,atomPos+diff);

% Phase due to propgation
if phase_delay==true
    delay_phase = zeros(size_array);
    for i=1:size_array(2)
        delay_phase(:,i) = neff*(i-1)*period/wavelength;
    end
elseif phase_delay==false
        phase_delay=0;
else
    error("Wrong input");
end


% Normalize the input phase data
Phase=NorPhase(Phase);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%      Set an breakpoint this line to check the phase data.               %
%                                                                         %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Truncate the phase data at desired interval
start_index=26;
stop_index=73;
Phase = Truncated_Phase(Phase,start_index,stop_index);
T = T(1,start_index:stop_index);
R_range = R_range(1,start_index:stop_index);

% Taking propagation phase into consideration
% Subtracting the phase due to propagation(equivalent to a plane wave under
% the meta-atoms)
tic; % timer start


% Create the desired phase profile
if lens_type=="axicon"
    Dphase = axiconOutput(delay_phase,phi0,beta,center,atomPos_X,atomPos_Y,1.55);
elseif lens_type=="spherical"
    Dphase = SphericalOutput(delay_phase,phi0,f,center,atomPos_X,atomPos_Y,1.55);
elseif lens_type=="grating"
    Dphase = gratingOutput(delay_phase,phi0,gt_angle,period,atomPos_X,atomPos_Y,1.55);
elseif lens_type=="None"
    Dphase = zeros(1,length(atomPos));
else
    error("Wrong lens type");
end
%Dphase2 = SphericalOutput(delay_phase2,phi0,f,[0,0],atomPos_X2,atomPos_Y2,1.55);

% Do the interpolation to find the corresponding radius and transmission data.
[R_list,T_list]=findSize(Dphase,Phase,T,R_range);
%[R_list2,T_list2]=findSize(Dphase2,Phase,T,R_range);


% Output radius list
if outputlist==true
    writematrix(R_list,outputname,'Delimiter','tab')
end

% Simulating energy decay below the waveguide
% base on exp(-alpha*x)
if decay==true
    for i=1:size_array(2)
        T_list(:,i) = T_list(:,i)*exp(log(decay_rate)/lens_size(1)*i*period);
    end
end
%T_list = T_list/max(max(T_list));
Phase_dist = Dphase+1;
% Check the Radius and Transmission distribution
figure;
surface(atomPos_X,atomPos_Y,T_list);
title("Top emission intensity distribution");

figure;
surface(atomPos_X,atomPos_Y,Phase_dist);
title("Phase distribution");


% Stop here if you only want to find the corresponding radii list of the
% designed phase profile.

%% Use superposition of the point source to estimate the result

x_range = input.x_range;
x_res = input.x_res;
y_range = input.y_range;
y_res = input.y_res;
z_range = input.z_range;
z_res = input.z_res;
% Calculating the field and plot it out.(real, imag, and abs)
if plot_focusingField==true
    focusing_field=Focusing_Slice(Dphase,T_list,atomPos,...
        x_range,z_range,0,x_res,z_res,wavelength,false);

    [a,b]=max(max(abs(focusing_field(:,100:end))));
    focal_z = z_range(1)+(b+100)*(z_range(2)-z_range(1))/z_res;
end


if plot_focalField==true
    focal_field=Focal_Slice(Dphase,T_list,atomPos,...
    x_range,y_range,focal_z,x_res,y_res,wavelength,false); 
end

% Find the maximum value from the focusing plot and show the focal slice of
% the plane.
%[a,b]=max(max(abs(focusing_field)));
%focal_z = z_range(1)+b*(z_range(2)-z_range(1))/z_res;
%focal_field=Focal_Slice(Dphase,T_list,atomPos,...
%    x_range,y_range,focal_z,x_res,y_res,1.55,false); 


toc; %timer stop








