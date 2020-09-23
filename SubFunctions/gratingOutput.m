function Dphase = gratingOutput(init_phase,phi0,angle,period,atomPos_X,atomPos_y,lambda)
% This function generates a phase list of meta-atoms that acts like an
% grating coupler.
% Unit in "micron".
% Phase profile: -(2pi/lambda)*(sqrt(r^2+f^2)-f)
% Input:
% init_phase = initial phase data (normalized data, and between[-1,0]!!)
%         Should be 2xN array.(phase of each atom on xy-plane)
%         If there is no initial phase data, just input '0'.
% phase = phase data, which is used to avoid empty region in the middle part.
% angle = grating angle (degree)
% mid_point = projection of the focal point on the xy-plane.
% atom_pos = position of each atom. Should be 2D array.(um)
% lambda = wavelength(um)
% Output
% list: radius list of spherical phase profile
% Dphase: 1xN corresponding phase array(Normalized, starting from -1)
if init_phase==0
    Dphase = zeros(1,length(atom_pos));
else
    Dphase = init_phase;
end
lens_size = size(atomPos_X);

for i=1:lens_size(2)
    Dphase(:,i) = Dphase(:,i)+i*period*sin(deg2rad(angle))/lambda+1+phi0;

end
for i=1:lens_size(1)
    for j=1:lens_size(2)
        while Dphase(i,j)<-1
            Dphase(i,j)=Dphase(i,j)+1;
        end
        while Dphase(i,j)>0
             Dphase(i,j)=Dphase(i,j)-1;
        end
    end
    
end

end

