function Dphase=axiconOutput(init_phase,phase,beta,mid_point,atomPos_X,atomPos_Y,lambda)
% This function generates a phase list of meta-atoms that acts like an axion.
% Unit in "micron".
% Phase profile: -(2pi/lambda)*(sqrt(r^2+f^2)-f)
% Input:
% init_phase = initial phase data (normalized data, and between[-1,0]!!)
%         Should be 2xN array.(phase of each atom on xy-plane)
%         If there is no initial phase data, just input '0'.
% phase = phase data, which is used to avoid empty region in the middle part.
% beta = beta of axicon (degree)
% mid_point = projection of the focal point on the xy-plane.
% atom_pos = position of each atom. Should be 2D array.(um)
% lambda = wavelength(um)
% Output
% list: radius list of spherical phase profile
% Dphase: corresponding phase array(Normalized, starting from -1)
% 
beta = beta*pi/180;
% Check if there are propagation angle need to be considered
if init_phase==0
    Dphase = zeros(size(atomPos_x));
else
    Dphase = init_phase;
end
lens_size = size(atomPos_X);
x0=mid_point(1);
y0=mid_point(2);
r_square = (atomPos_X-x0).^2+(atomPos_Y-y0).^2;
Dphase = Dphase-sqrt(r_square)*sin(beta)/lambda;
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