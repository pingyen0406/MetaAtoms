function Dphase=SphericalOutput(init_phase,phi0,f,mid_point,atomPos_X,atomPos_Y,lambda)
% This function generates a phase list of meta-atoms that can focus the beam.
% Unit in "micron".
% Phase profile: -(2pi/lambda)*(sqrt(r^2+f^2)-f)
% Input:
% init_phase = initial phase data (normalized data, and between[-1,0]!!)
%         Should be 2xN array.(phase of each atom on xy-plane)
%         If there is no initial phase data, just input '0'.
% phase = phase data, which is used to avoid empty region in the middle part.
% f = focal length(um)
% mid_point = projection of the focal point on the xy-plane.
% atom_pos = position of each atom. Should be 2D array.(um)
% lambda = wavelength(um)
% Output
% list: radius list of spherical phase profile
% Dphase: 1xN corresponding phase array(Normalized, starting from -1)
% 

% Calculate the desired phase profile.
if init_phase==0
    Dphase = zeros(size(atomPos_X));
else
    Dphase = init_phase;
end
x=mid_point(1);
y=mid_point(2);
lens_size = size(atomPos_X);
r_square = (atomPos_X-x).^2+(atomPos_Y-y).^2;
Dphase = -Dphase-(sqrt(r_square+f^2)-f)/lambda+phi0;
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