function Dphase=SphericalOutput(init_phase,f,mid_point,atom_pos,lambda)
% This function generates a phase list of meta-atoms that can focus the beam.
% Unit in "micron".
% Phase profile: -(2pi/lambda)*(sqrt(r^2+f^2)-f)
% Input:
% init_phase = initial phase data (normalized data, and between[-1,0]!!)
%         Should be 2xN array.(phase of each atom on xy-plane)
%         If there is no initial phase data, just input '0'.
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
    Dphase = zeros(1,length(atom_pos));
else
    Dphase = init_phase;
end
x=mid_point(1);
y=mid_point(2);
for i=1:length(atom_pos)
    r_square = (atom_pos(1,i)-x)^2+(atom_pos(2,i)-y)^2;
    Dphase(1,i) = Dphase(1,i)-(sqrt(r_square+f^2)-f)/lambda;
    while Dphase(1,i)<-1
        Dphase(1,i)=Dphase(1,i)+1;
    end
end
end