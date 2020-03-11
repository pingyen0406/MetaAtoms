function field=Focusing_Slice(Phase,T,atom_pos,x_range,z_range,y,x_res,z_res,lambda)
% Calaulating and plotting the focusing behavior from the meta-atoms to the
% focal spot.
% Unit in "micron".
% E = (1/r)*A*exp(i*k*r)*exp(i*starting phase)
% Phase: starting phase at z=0, and should be 1D array.
% atom_pos: x,y position array of meta-atoms, which should be 2xN array.
% T : transmission = amplitude, and should be 1D array.
% x_range & z_range = [xmin,xmax]
% x_res & z_res: resolution in x & y.
% y: where to slice
% lambda: operating wavelength 
x_list = linspace(x_range(1),x_range(2),x_res);
z_list = linspace(z_range(1),z_range(2),z_res);
r_list = zeros(length(x_list),length(z_list));
field = zeros(length(x_list),length(z_list),length(T));
% r_list is a 3-dim matrix. r_list(i,j,k) means distance from k-th 
% meta-atom to x_list(i) and z_list(j). So is field.
for i= 1:length(x_list)
    for j=1:length(z_list)
        for k=1:length(atom_pos)
            tmp=(x_list(i)-atom_pos(1,k))^2+(y-atom_pos(2,k))^2+(z_list(j))^2;
            r_list(i,j,k)=sqrt(tmp);
        end
    end
end
prop_phase = r_list/lambda;
for i=1:length(x_list)
    for j=1:length(z_list)
        for k=1:length(T)
            field(i,j,k) = (1/r_list(i,j,k))*T(k)*...
                exp(1i*2*pi*prop_phase(i,j,k))*exp(1i*2*pi*Phase(k));
        end
    end
end
field = sum(field,3);
figure(1);
colormap('jet');
image(z_list,x_list,real(field),'CDataMapping','scaled');
title(['real part at y=',num2str(y),'um']);
figure(2);
colormap('jet');
image(z_list,x_list,imag(field),'CDataMapping','scaled');
title(['imag part at y=',num2str(y),'um']);
figure(3);
colormap('jet');
image(z_list,x_list,abs(field),'CDataMapping','scaled');
title(['Absolute value at y=',num2str(y),'um']);
end