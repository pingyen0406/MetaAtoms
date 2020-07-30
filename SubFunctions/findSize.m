function [R_list, T_list] = findSize(Dphase,Phase,T,R)
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
lensSize = size(Dphase);
R_list = zeros(size(Dphase));
for i =1:lensSize(1)
    for j=1:lensSize(2)
        if isnan(interp1(Phase,R,Dphase(i,j)))==true
            R_list(i,j)=0;
        else
            R_list(i,j) = interp1(Phase,R,Dphase(i,j));
        end
    end
end
T_list = zeros(size(Dphase));
for i =1:lensSize(1)
    for j=1:lensSize(2)
        if isnan(interp1(R,T,R_list(i,j)))==true
            T_list(i,j)=0;
        else
            T_list(i,j) = interp1(R,T,R_list(i,j));
        end
    end
end
end