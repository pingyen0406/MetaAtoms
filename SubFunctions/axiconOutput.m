function Dphase=axiconOutput(init_phase,beta,mid_point,atom_pos,lambda)
% This function generates a phase list of meta-atoms that acts like an axion.
% Unit in "micron".
% Phase profile: -(2pi/lambda)*(sqrt(r^2+f^2)-f)
% Input:
% init_phase = initial phase data (normalized data, and between[-1,0]!!)
%         Should be 2xN array.(phase of each atom on xy-plane)
%         If there is no initial phase data, just input '0'.
% beta = beta of axicon (degree)
% mid_point = projection of the focal point on the xy-plane.
% atom_pos = position of each atom. Should be 2D array.(um)
% lambda = wavelength(um)
% Output
% list: radius list of spherical phase profile
% Dphase: 1xN corresponding phase array(Normalized, starting from -1)
% 
beta = beta*pi/180;
% Check if there are propagation angle need to be considered
if init_phase==0
    Dphase = zeros(1,length(atom_pos));
else
    Dphase = init_phase;
end
x0=mid_point(1);
y0=mid_point(2);
for i=1:length(atom_pos)
    r_square = (atom_pos(1,i)-x0)^2+(atom_pos(2,i)-y0)^2;
    Dphase(1,i) = Dphase(1,i)-sqrt(r_square)*sin(beta)/lambda;
    while Dphase(1,i)<-1
        Dphase(1,i)=Dphase(1,i)+1;
    end
end
end