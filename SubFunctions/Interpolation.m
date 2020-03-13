function [R_list, T_list] = Interpolation(Dphase,Phase,T,R)
% Find the corresponding radius and transmission of the meta-atoms by 
% interpolation.
% Find the radius list first and then use the radius list to do
% interpolation in order to have the corresponding transmission list.
% Input:
% Dphase: designed phase profile
% Phase: phase profile of all the radius list
% T: transmission value of all the radius list
% R: radius list
% Output:
% R_list: radius list comes from designed phase profile
% T_list: transmission list comes from designed radius list.
N = length(Dphase);
R_list = zeros(1,N);
for i =1:N
    if isnan(interp1(Phase,R,Dphase(i)))==1
        R_list(1,i)=0;
    else
        R_list(1,i) = interp1(Phase,R,Dphase(i));
    end
end
T_list = zeros(1,N);
for i=1:N
    if isnan(interp1(R,T,R_list(1,i)))==1
        T_list(1,i)=0;
    else
        T_list(1,i) = interp1(R,T,R_list(1,i));
    end
end
end