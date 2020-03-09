function [list,Dphase]=SphericalOutput(Phase,T,R,f,lattice,N,lambda)
% This function generates a list of meta-atoms that can focus the beam.
% Note: The initial phase value of every points is considered identically.
% Unit in "micron".
% Phase profile: -(2pi/lambda)*(sqrt(r^2+f^2)-f)
% Input:
% Phase = phase data (Normalized data, i.e. between[-1,1])
% T = transmission data
% R = radius list(um)
% f = focal length(um)
% lattice = lattice constant(um)
% N = number of atoms
% lambda = wavelength(um)
% Output
% list: radius list of spherical phase profile
% Dphase: corresponding phase of list(Normalized)

% Calculate the desired phase profile.
Dphase = zeros(1,N);
list = zeros(1,N);
for i=1:N
    Dphase(1,i) = -(sqrt((abs(i-1-floor(N/2))*lattice)^2+f^2)-f)/lambda;
    while Dphase(1,i)<-1
        Dphase(1,i)=Dphase(1,i)+1;
    end
end
% Find the corresponding radius of the meta-atoms by interp.
for i =1:N
    if isnan(interp1(Phase,R,Dphase(i)))==1
        list(1,i)=0;
    else
        list(1,i) = interp1(Phase,R,Dphase(i));
    end
end
end