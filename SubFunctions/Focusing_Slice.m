function field=Focusing_Slice(Phase,T,atomPos,x_range,z_range,y,x_res,z_res,lambda,symmetry)
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
% symmetry: true or false. Use symmetric property to reduce time.
type = "focusing";
x_list = linspace(x_range(1),x_range(2),x_res);
z_list = linspace(z_range(1),z_range(2),z_res);

if symmetry == false
    field = pointSource_field(Phase,T,atomPos,x_list,z_list,y,1.55,type);
elseif symmetry == true
    x_len = floor(length(x_list)/2);
    field = pointSource_field(Phase,T,atomPos,x_list(1:x_len),...
        z_list,y,lambda,type);
    field = cat(1,field,flip(field,1));
    % If the interval is odd, it needs to calculate the middle line.
    if mod(x_len,2)~=0
        mid_field = pointSource_field(Phase,T,atomPos,x_list(x_len+1),...
        z_list,y,lambda,type);
        field = [field(1:x_len,:);mid_field;field(x_len+1:end,:)];
    end
else
    error("Wrong input on symmetry");
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