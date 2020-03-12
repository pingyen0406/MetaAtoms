function phi = PropCorrect(N,period,neff,wavelength)
% Correcting the desire phase list by taking propagation phase into
% consideration.
% Input:
% period: period of the metasurface
% neff: effective index while propagating through the waveguide
% wavelength: operating wavelength
% Output:
% phi: phase profile of the propagating wave(normalized,between[0,1])
phi = zeros(1,N);
for i=1:N
    phi(1,i)=-(i-1)*period*neff/wavelength;    
    while phi(1,i)<-1
        phi(1,i)=phi(1,i)+1;
    end
end

end