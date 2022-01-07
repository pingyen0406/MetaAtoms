% Calculating the scalar field emitted by 2D meta-atom array. 
% General using method:
% 1. Setting parameters in the input.json(transmission & phase data, lens
% size, radius & height range and step)
% 2. For the 1st run, set a breakpoint at line 60. You need to choose the
% range of the radii.
% 3. After determining the radii range, just follow the pop out dialog.
% Note that if you choose "custom", the input desired phase matrix should
% have value between 0 to 2pi.



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
rowIndex = (height-input.heightRange(1))/input.heightStep+1;
T = all_T(rowIndex,:);
Phase = all_Phase(rowIndex,:);
R_range = [input.radiusRange(1):input.radiusStep:input.radiusRange(2)];

% Parameters
period = input.period;
lens_size = input.size;
% 'Radius' or 'Length' of metalens (circle or square)
center=input.center; % middle point of the pattern
neff = input.neff; % effective index derived from FDTD
wavelength = 1.55;
use_prop_phase=input.prop_phase;
decay = input.dacay;
decay_rate = 0.5; % The ramaining power in the waveguide. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[atomPos_X,atomPos_Y] = squarePos([0,0],period,lens_size);
size_array = size(atomPos_X);

%diff = [ones(1,N*N)*N*period;zeros(1,N*N)];
% Use 3 squares to create a rectangle(20um*60um lens)
%atomPos = cat(2,atomPos-diff,atomPos,atomPos+diff);

% Normalize the input phase data
Phase=NorPhase(Phase,true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%      Set an breakpoint this line to check the phase data.               %
%                                                                         %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot(R_range,Phase);
% Truncate the phase data at desired interval
start_index=26;
stop_index=73;
Phase = Truncated_Phase(Phase,start_index,stop_index);
T = T(1,start_index:stop_index);
R_range = R_range(1,start_index:stop_index);


% Taking propagation phase into consideration
% Phase due to propgation in the wg
if use_prop_phase==true
    prop_phase = zeros(size_array);
    for i=1:size_array(2)
        prop_phase(:,i) = neff*(i-1)*period/wavelength;
    end
elseif use_prop_phase==false
       prop_phase=0;
else
    error("Wrong input");
end

tic; % timer start

% Pop out dialog to choose lens type
lensTypes = {'axicon','spherical','custom','none'};
[indx,~] = listdlg('PromptString','Select a type of lens.',...
    'SelectionMode','single','ListString',lensTypes);
lens_type=lensTypes(indx);

% Create the desired phase profile
if lens_type=="axicon"
    prompt = {'Enter beta(degree):'};
    dlgtitle = 'Input';
    dims = 1;
    definput = {'5'};
    beta = inputdlg(prompt,dlgtitle,dims,definput);
    beta = str2double(beta{1});
    phi0 = askPhi0;
    Dphase = axiconOutput(prop_phase,phi0,beta,center,atomPos_X,atomPos_Y,1.55);
elseif lens_type=="spherical"
    prompt = {'Enter focal length(um):'};
    dlgtitle = 'Input';
    dims = 1;
    definput = {'100'};
    f = inputdlg(prompt,dlgtitle,dims,definput);
    f = str2double(f{1});
    phi0 = askPhi0;
    Dphase = SphericalOutput(prop_phase,phi0,f,center,atomPos_X,atomPos_Y,1.55);
elseif lens_type=="grating"
    prompt = {'Enter grating angle(degree):'};
    dlgtitle = 'Input';
    dims = 1;
    definput = {'0'};
    gt_angle = inputdlg(prompt,dlgtitle,dims,definput);
    gt_angle = str2double(gt_angle{1});
    phi0 = askPhi0;
    Dphase = gratingOutput(prop_phase,phi0,gt_angle,period,atomPos_X,atomPos_Y,1.55);
elseif lens_type=="custom"
    customText = ['Choose the custom input file. Matrix size = ', num2str(size_array(1)),'x',...
        num2str(size_array(2)),'.'];
    disp(customText);
    [customFile,customPath] = uigetfile('*.txt');
    hologram = NorPhase(readmatrix([customPath,customFile]),true);
    phi0 = askPhi0;
    Dphase = customOutput(prop_phase,hologram,phi0,atomPos_X,atomPos_Y);
elseif lens_type=="None"
    Dphase = zeros(size_array(1),size_array(2));
else
    error("Wrong lens type");
end



% Do the interpolation to find the corresponding radius and transmission data.
[R_list,T_list,Phase_list]=findSize(Dphase,Phase,T,R_range);
%writematrix(Phase_list,"1D_gt_n20_phase.txt");
% Calculate the average radius
all_r=0;
for i=1:size_array(1)
    for j=1:size_array(2)
        all_r = all_r+R_list(i,j);
    end
end
avg_r = all_r/size_array(1)/size_array(2);

for i=1:size_array(1)
    for j=1:size_array(2)
        if j==1
            if R_list(i,j)==0
                R_list(i,j) = R_list(i-1,j);
            end
        else
            if R_list(i,j)==0
                R_list(i,j) = R_list(i,j-1);
            end
        end
    end
end
% Output radius list
if outputlist==true
    writematrix(R_list,outputname,'Delimiter','space');
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
colormap('jet');
figure;
surface(atomPos_X,atomPos_Y,R_list);
title("Radius distribution");


% Function of asking phi0

%% Use superposition of the point source to estimate the result
% Not using for a long time, it may need some modifacation.
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
function phi0 = askPhi0()
prompt = {'Enter phi0 phase shift:'};
dlgtitle = 'Input';
dims = 1;
definput = {'0'};
phi0 = inputdlg(prompt,dlgtitle,dims,definput);
phi0 = str2double(phi0{1});
end







