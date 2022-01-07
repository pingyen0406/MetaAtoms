function Dphase = customOutput(init_phase,customPhase,phi0,atomPos_X,atomPos_Y)
% This function reads an input custom phase and calibrates the propagation phase.
% Unit in "micron".
% Input:
% init_phase = initial propagation phase (normalized data, and between[-1,0]!!)
%         Should be 2xN array.(phase of each atom on xy-plane)
%         If there is no initial phase data, just input '0'.
% customPhase = customPhase data between[-1,0]!!
% % atomPos_X & atomPos_Y = position of each atom. Should be 2D array.(meshgrid)
% lambda = wavelength(um)
% Output
% Dphase: Phase matrix on the x-y grids(Normalized, starting from -1)
% 
if size(customPhase)~=size(atomPos_X)
    error('customPhase input size does not equal to metalens size!')
end
if init_phase==0
    Dphase = zeros(size(atomPos_X));
else
    Dphase = -init_phase;
end
% Calculate the desired phase profile.
lens_size = size(atomPos_X);
Dphase = Dphase+customPhase+phi0;

% Wrapping the phase data
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


