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
fname_T = [folder_path,'SweepT562.txt'];
fname_Phase = [folder_path,'SweepPhase562.txt'];




plot_focalField = input.plot_focalField;
plot_focusingField = input.plot_focusingField;
outputlist = input.outputlist;
outputname = [folder_path,'Al2O3test_2.txt'];

% Read data from S4 calculation
% (i,j) is 500+10*i nm height and 100+1*j nm radius
all_T = readmatrix(fname_T);
all_Phase = readmatrix(fname_Phase);
height = 1200;
T = all_T(height/10-50+1,:);
Phase = all_Phase(height/10-50+1,:);
R_range = [0.03:0.002:0.248];

% Parameters
period = input.period;
lens_type = input.lens_type; % "axicon" or "spherical"
size = input.size;
% 'Radius' or 'Length' of metalens (circle or square)
f_num = 3.5; % f-number
f = input.f; % focal length
center=input.center; % middle point of the pattern
beta =5; % beta angle of axicon(in degree)
neff = 2.858; % effective index derived from FDTD
wavelength = 1.55;
phase_delay=input.phase_delay;
decay = input.dacay;
decay_rate = 0.5;
x_range = input.x_range;
x_res = input.x_res;
y_range = input.y_range;
y_res = input.y_res;
z_range = input.z_range;
z_res = input.z_res;

%atomPos = squarePos("circle",[0,0],period,size/2);
atomPos = squarePos("square",[0,0],period,size);
N = sqrt(length(atomPos));
diff = [ones(1,N*N)*N*period;zeros(1,N*N)];
% Use 3 squares to create a rectangle
atomPos = cat(2,atomPos-diff,atomPos,atomPos+diff);

%---------------------------- Test data-----------------------------------
%{
tic;
test_name = 'perfectAtom.txt';
test_f = readmatrix(test_name);
test_r = transpose(test_f(:,1));
test_T = transpose(test_f(:,2));
test_Phase = NorPhase(transpose(test_f(:,3)));
test_Phase(1,end)=0;
Dphase = SphericalOutput(0,test_Phase,f,center,atomPos,1.55);
[R_list,T_list]=Interpolation(Dphase,test_Phase,test_T,test_r);
focusing_field=Focusing_Slice(Dphase,T_list,atomPos,...
        x_range,z_range,0,x_res,z_res,1.55,true);
focal_field=Focal_Slice(Dphase,T_list,atomPos,...
    x_range,y_range,f,x_res,y_res,1.55,true); 
toc;
%}
%--------------------------------------------------------------------------



% Normalize the input phase data
Phase=NorPhase(Phase);
% Set an breakpoint this line to check the phase data.
% Truncate the phase data at desired interval
start_index=35;
stop_index=79;
Phase = Truncated_Phase(Phase,start_index,stop_index);
T = T(1,start_index:stop_index);
R_range = R_range(1,start_index:stop_index);

% Taking propagation phase into consideration
% Subtracting the phase due to propagation(equivalent to a plane wave under
% the meta-atoms)
if phase_delay==true
    count=1;
    delay_phase = zeros(1,length(atomPos));
    prop_N = length(atomPos)/N;
    tmp_delay_phase = PropCorrect(prop_N,period,neff,wavelength); 
    tmp_delay_phase = NorPhase(tmp_delay_phase);
    for i=1:prop_N
        for j=1:N
            delay_phase(1,count)=tmp_delay_phase(1,i);
            count=count+1;
        end
    end
else
    delay_phase=0;
end


tic; % timer start

% Output multiple R_list at a time
%{
count=1;
for center=[[-30;0],[-10;0],[10;0],[30;0]]
    outputname = [folder_path,'focus_180_',num2str(count),'.txt'];
    Dphase = SphericalOutput(delay_phase,Phase,f,center,atomPos,1.55);
% Do interpolation to find the corresponding radius and transmission data.
    [R_list,T_list]=Interpolation(Dphase,Phase,T,R_range);
    
    outf = fopen(outputname,'w');
    for i=1:length(R_list)
        fprintf(outf,'%f\n',R_list(i));
    end
    fclose(outf);
    count=count+1;
end
%}
    
% Create the desired phase profile
if lens_type=="axicon"
    Dphase = axiconOutput(delay_phase,Phase,beta,center,atomPos,1.55);
elseif lens_type=="spherical"
    Dphase = SphericalOutput(delay_phase,Phase,f,center,atomPos,1.55);
else
    error("Wrong lens type");
end
% Do the interpolation to find the corresponding radius and transmission data.
[R_list,T_list]=Interpolation(Dphase,Phase,T,R_range);



% Output radius list
if outputlist==true 
    outf = fopen(outputname,'w');
    for i=1:length(R_list)
        fprintf(outf,'%f\n',R_list(i));
    end
    fclose(outf);
end

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
%figure;
%scatter3(atomPos(1,:),atomPos(2,:),T_list,'filled');
%title("Transmission distribution");
%figure;
%scatter3(atomPos(1,:),atomPos(2,:),Dphase,'filled');
%title("Phase distribution");


% Calculating the field and plot it out.(real, imag, and abs)
if plot_focusingField==true
   focusing_field=Focusing_Slice(Dphase,T_list,atomPos,...
        x_range,z_range,0,x_res,z_res,wavelength,false);
    
end
[a,b]=max(max(abs(focusing_field(:,100:end))));
focal_z = z_range(1)+(b+100)*(z_range(2)-z_range(1))/z_res;

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







