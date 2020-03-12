function field=Focal_Slice(Phase,T,atom_pos,x_range,y_range,z,x_res,y_res,lambda)
% Calaulating and plotting the focusing behavior on the focal plane.
% Unit in "micron".
% E = (1/r)*A*exp(i*k*r)*exp(i*starting phase)
% Phase: starting phase at z=0, and should be 1D array.
% atom_pos: x,y position array of meta-atoms, which should be 2xN array.
% T : transmission = amplitude, and should be 1D array.
% x_range & y_range = [xmin,xmax]
% x_res & y_res: resolution in x & y.
% z: where to slice
% lambda: operating wavelength 
x_list = linspace(x_range(1),x_range(2),x_res);
y_list = linspace(y_range(1),y_range(2),y_res);
r_list = zeros(length(x_list),length(y_list));
field = zeros(length(x_list),length(y_list),length(T));
% r_list is a 3-dim matrix. r_list(i,j,k) means distance from k-th 
% meta-atom to x_list(i) and y_list(j). So is field.
for i= 1:length(x_list)
    for j=1:length(y_list)
        for k=1:length(atom_pos)
            tmp=(x_list(i)-atom_pos(1,k))^2+(y_list(1,j)-atom_pos(2,k))^2+z^2;
            r_list(i,j,k)=sqrt(tmp);
        end
    end
end
prop_phase = r_list/lambda;
for i=1:length(x_list)
    for j=1:length(y_list)
        for k=1:length(T)
            field(i,j,k) = (1/r_list(i,j,k))*T(k)*...
                exp(1i*2*pi*prop_phase(i,j,k))*exp(1i*2*pi*Phase(k));
        end
    end
end
field = sum(field,3);
figure(1);
colormap('jet');
image(y_list,x_list,real(field),'CDataMapping','scaled');
title(['real part at z=',num2str(z),'um']);
figure(2);
colormap('jet');
image(y_list,x_list,imag(field),'CDataMapping','scaled');
title(['imag part at z=',num2str(z),'um']);
figure(3);
colormap('jet');
image(y_list,x_list,abs(field),'CDataMapping','scaled');
title(['Absolute value at z=',num2str(z),'um']);
end
