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
field = zeros(length(x_list),length(y_list));
tmpField=zeros(1,length(atom_pos));

for i= 1:length(y_list)
    for j=1:length(x_list)
        for k=1:length(atom_pos)
            tmp_rr=(x_list(1,j)-atom_pos(1,k))^2+(y_list(1,i)...
                -atom_pos(2,k))^2+z^2;
            tmp_phase=sqrt(tmp_rr)/lambda;    
            tmpField(k) = (1/sqrt(tmp_rr))*T(k)*...
                exp(1i*2*pi*tmp_phase)*exp(1i*2*pi*Phase(k));
        end
        field(i,j)=sum(tmpField);
    end
end
figure;
colormap('jet');
image(x_list,y_list,real(field),'CDataMapping','scaled');
title(['real part at z=',num2str(z),'um']);
xlabel('x');ylabel('y');
figure;
colormap('jet');
image(x_list,y_list,imag(field),'CDataMapping','scaled');
title(['imag part at z=',num2str(z),'um']);
xlabel('x');ylabel('y');
figure;
colormap('jet');
image(x_list,y_list,abs(field),'CDataMapping','scaled');
title(['Absolute value at z=',num2str(z),'um']);
xlabel('x');ylabel('y');
end
