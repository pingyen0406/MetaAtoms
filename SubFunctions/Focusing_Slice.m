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
field = zeros(length(x_list),length(z_list));
tmpField=zeros(1,length(atom_pos));

for i= 1:length(x_list)
    for j=1:length(z_list)
        for k=1:length(atom_pos)
            tmp_rr=(x_list(i)-atom_pos(1,k))^2+(y-atom_pos(2,k))^2+z_list(j)^2;
            tmp_phase=sqrt(tmp_rr)/lambda;    
            tmpField(k) = (1/sqrt(tmp_rr))*T(k)*...
                exp(1i*2*pi*tmp_phase)*exp(1i*2*pi*Phase(k));
        end
        field(i,j)=sum(tmpField);
    end
end
figure;
colormap('jet');
image(z_list,x_list,real(field),'CDataMapping','scaled');
title(['real part at y=',num2str(y),'um']);
xlabel('z');ylabel('x');
figure;
colormap('jet');
image(z_list,x_list,imag(field),'CDataMapping','scaled');
title(['imag part at y=',num2str(y),'um']);
xlabel('z');ylabel('x');
figure;
colormap('jet');
image(z_list,x_list,abs(field),'CDataMapping','scaled');
title(['Absolute value at y=',num2str(y),'um']);
xlabel('z');ylabel('x');
end