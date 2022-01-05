clear all;

inf_P = "Lib562/80nmSiO2/newNeff/range2/SweepPhase562.txt";
inf_T = "Lib562/80nmSiO2/newNeff/range2/SweepT562.txt";
inf_Rlist = "Grating/80nmSiO2_1D_gt_20.txt";
all_Phase = readmatrix(inf_P);
all_T = readmatrix(inf_P);
R_list = readmatrix(inf_Rlist);
R_list=R_list(1,:);
angle=20;
neff=2.85;
period=0.562;
wavelength=1.55;
size_array = size(R_list);
delay_phase = zeros(size_array);
for i=1:size_array(2)
    delay_phase(i) = neff*(i-1)*period/wavelength;
end
%% Retrace the phase(from radii to phase)
height = 1.2;
Phase_lib = all_Phase(int16(height/0.05)-8+1,:);
R_range = [0.1:0.002:0.194];
Phase_lib = Phase_lib(26:73)/2/pi;
Phase_lib=NorPhase(Phase_lib,false);
lensSize = size(R_list);
Phase_list = zeros(size(R_list));

for i =1:lensSize(1)
    for j=1:lensSize(2)
        
        if isnan(interp1(R_range,Phase_lib,R_list(i,j)))==true
            Phase_list(i,j)=-1;
        else
            Phase_list(i,j) = interp1(R_range,Phase_lib,R_list(i,j));
        end
    end
end
writematrix(Phase_list,"1D_gt_20.txt")

%delay_phase = NorPhase(delay_phase);

%% Compare different height
for height=[1.2]

    T = all_T(int16(height/0.05)-8+1,:);
    Phase_lib = all_Phase(int16(height/0.05)-8+1,:);
    R_range = [0.1:0.002:0.194];
    Phase_lib = Phase_lib(26:73)/2/pi;
    Phase_lib=NorPhase(Phase_lib,false);
    lensSize = size(R_list);
    Phase_list = zeros(size(R_list));
    idealPhase = zeros(size(R_list)); 
    count=0;
    for i =1:lensSize(1)
        for j=1:lensSize(2)
            idealPhase(i,j) = (j-1)*period*sin(deg2rad(angle))/wavelength;
            if isnan(interp1(R_range,Phase_lib,R_list(i,j)))==true
                count=count-1;
                
                Phase_list(i,j)=-1;
            else
                Phase_list(i,j) = interp1(R_range,Phase_lib,R_list(i,j));
            end

        end
    end
    atom_phase = Phase_list;
    figure(1);plot(atom_phase,'LineWidth',1.5);
    hold on
    for i=1:lensSize(1)
        for j=1:lensSize(2)
            if Phase_list(i,j)~=-1
                 Phase_list(i,j) = Phase_list(i,j)+delay_phase(i,j);
%                 while Phase_list(i,j)>=0
%                     Phase_list(i,j)=Phase_list(i,j)-1;
%                 end
%                 while Phase_list(i,j)<-1
%                     Phase_list(i,j)=Phase_list(i,j)+1;
%                 end
           end
       end
   end
    %Phase_list=NorPhase(Phase_list);
    figure(2);
    plot(Phase_list(1:30),'LineWidth',1.5);
    
    hold on
end
% figure(1);
% hold on 
% plot(idealPhase,'LineWidth',1.5);
% figure(2);
% hold on
% plot((NorPhase(delay_phase)),'LineWidth',1.5);
figure(3);
plot(NorPhase(delay_phase,false));
hold on 
plot(NorPhase(idealPhase,false));
plot(NorPhase(idealPhase-delay_phase,false));

