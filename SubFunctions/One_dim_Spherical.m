function Dphase=One_dim_Spherical(Phase,f,lattice,N,lambda)
% This function generates a 1D phase profile of meta-atoms that can focus the beam.
% % Unit in "micron".
% Phase profile: -(2pi/lambda)*(sqrt(r^2+f^2)-f)
% Input:
% Phase = initial phase data (normalized data, and between[-1,0]!!)
%         If there is no initial phase data, just input '0'.
% f = focal length(um)
% lattice = lattice constant(um)
% N = number of atoms
% lambda = wavelength(um)
% Output
% list: radius list of spherical phase profile
% Dphase: corresponding phase of list(Normalized, starting from -1)

% Calculate the desired phase profile.
if Phase==0
    Dphase = Phase;
else
    Dphase = zeros(1,N);
end

for i=1:N
    Dphase(1,i) = Dphase(1,i)-(sqrt((abs(i-1-floor(N/2))*lattice)^2+f^2)-f)/lambda;
    while Dphase(1,i)<-1
        Dphase(1,i)=Dphase(1,i)+1;
    end
end

end