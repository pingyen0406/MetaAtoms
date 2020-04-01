function field=Focal_Slice(Phase,T,atomPos,x_range,y_range,z,x_res,y_res,...
    lambda,symmetry)
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
% symmetry: true or flase. Use x-y symmetric property to reduce 3/4 time.
x_list = linspace(x_range(1),x_range(2),x_res);
y_list = linspace(y_range(1),y_range(2),y_res);
type = "focal";

if symmetry==false
    field = pointSource_field(Phase,T,atomPos,x_list,y_list,z,lambda,type);
      
elseif symmetry==true  
    x_len = floor(length(x_list)/2);
    y_len = floor(length(y_list)/2);
    field = pointSource_field(Phase,T,atomPos,x_list(1:x_len),...
        y_list(1:y_len),z,lambda,type);
    field = cat(2,field,flip(field,2));
    field = cat(1,field,flip(field,1));
    % If the interval is odd, it needs to calculate the middle line.
    if mod(length(x_list),2)~=0
        mid_field1=pointSource_field(Phase,T,atomPos,x_list(x_len+1),...
            y_list(1:y_len),z,lambda,type);
        center_field = pointSource_field(Phase,T,atomPos,x_list(x_len+1),...
            y_list(y_len+1),z,lambda,type);
        mid_field_1 = cat(1,midfield_1,flip(mid_field1,1)); % vertical field without center point
        mid_field_2 = [mid_field_1(1,1:x_len)',center_field,mid_field_1(1,x_len+1:end)'];
        field = [field(:,1:y_len),mid_field_1,field(:,y_len+1:end)];
        field = [field(1:x_len,:);mid_field_2;field(x_len+1:end,:)];
    end
     
else
    error("Wrong input on symmetry");
end
field = field';

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
