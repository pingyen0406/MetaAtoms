function field_matrix = pointSource_field(Phase,T,atomPos_X,atomPos_Y,interval_1,...
    interval_2,slicePoint,lambda,type)
% Calaulating the 2D field slice which emitted from meta-atoms at given position.
% Unit in "micron".
% Input:
% E = (1/r)*A*exp(i*k*r)*exp(i*starting phase)
% Phase: starting phase at z=0, and should be 1D array.
% atom_pos: x,y position array of meta-atoms, which should be 2xN array.
% T : transmission = amplitude, and should be 1D array.
% interval_1: grid point array of x (plane of metasurfaces)
% interval_2: grid point array of y or z
% slicePoint: where to slice on normal axis
% lambda: operating wavelength 
% Output:
% field_matrix: scalar value on given intervals
% type: 'focal' or 'focusing' plane. The distance term has slight
% difference.
tmpField=zeros(size(atomPos_X));
field_matrix = zeros(length(interval_1),length(interval_2));
if type == "focal"
    for i= 1:length(interval_1)
        for j=1:length(interval_2)
            for k=1:length(atomPos)
                tmp_rr=(interval_1(1,i)-atomPos(1,k))^2+(interval_2(1,j)...
                    -atomPos(2,k))^2+slicePoint^2;
                tmp_phase=sqrt(tmp_rr)/lambda;    
                tmpField(k) = (1/sqrt(tmp_rr))*T(k)*...
                    exp(1i*2*pi*tmp_phase)*exp(1i*2*pi*Phase(k));
            end
            field_matrix(i,j)=sum(tmpField);
        end
    end
elseif type== "focusing"
    for i= 1:length(interval_1)
        for j=1:length(interval_2)
            for k=1:length(atomPos)
                tmp_rr=(interval_1(1,i)-atomPos(1,k))^2+(slicePoint...
                    -atomPos(2,k))^2+interval_2(1,j)^2;
                tmp_phase=sqrt(tmp_rr)/lambda;    
                tmpField(k) = (1/sqrt(tmp_rr))*T(k)*...
                    exp(1i*2*pi*tmp_phase)*exp(1i*2*pi*Phase(k));
            end
            field_matrix(i,j)=sum(tmpField);
        end
    end
else
    error('Wrong input of type option');
end
end

