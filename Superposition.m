clear all;clc;
% For Linux path
%folder_path = '/home/pingyen/Simulation/MATLAB/MetaAtoms/Lib562/60nmAl2O3/';
%addpath("/home/pingyen/Simulation/MATLAB/Lib562/MetaAtoms/");
% For Windows path
folder_path = 'D:/Dropbox/MATLAB/MetaAtoms/Lib562/60nmAl2O3/TopAl2O3/';
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
lattice = 0.562;
f = 46.8; % focal length
N = 100; % number of meta-atoms
neff = 2.858; % effective index derived from FDTD
atomPos = zeros(2,N);
% crate x-position array of meta-atoms
for i=1:N
    atomPos(1,i)=lattice*i;
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
% Choosing desired part of phase data

Phase=NorPhase(Phase);
%plot(R_list,Phase);
% Set an breakpoint this line to check the phase data.

Phase = Truncated_Phase(Phase,32,79);
T = T(1,32:79);
R_list = R_list(1,32:79);
Dphase = PropCorrect(period,neff,wavelength);
Dphase = NorPhase(Dphase);
Dphase = SphericalOutput(Dphase,f,lattice,N,1.55);
[R_list,T_list]=Interpolation(Dpahse,T,R_list);
plot(T_list);
hold on
if plot_field==true
    Field=Eatom(Dphase,T_list,atomPos,[0,60],[1,61],0,600,600,1.55);
end

% Output List
if outputlist==true 
    outf = fopen(outputname,'w');
    for i=1:length(sphericalList)
        fprintf(outf,'%f\n',sphericalList(i));
    end
    fclose(outf);
end


