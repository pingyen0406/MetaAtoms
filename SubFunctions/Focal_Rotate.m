function field=Focal_Rotate(Phase,T,atomPos,x_range,y_range,z,x_res,y_res,...
    lambda)
% Calaulating and plotting the focusing behavior on the focal plane. Using
% rotational symmetry.
% Unit in "micron".
% E = (1/r)*A*exp(i*k*r)*exp(i*starting phase)
% Phase: starting phase at z=0, and should be 1D array.
% atom_pos: x,y position array of meta-atoms, which should be 2xN array.
% T : transmission = amplitude, and should be 1D array.
% x_range & y_range = [xmin,xmax]
% x_res & y_res: resolution in x & y.
% z: where to slice
% lambda: operating wavelength 
field=zeros(0,0);
type = "focal";
      
% First, calculating a single line from x,y=0,0 to x,y=x_max,0.
% And then transforming it into polar system and rotating it.
% After thet, transforming it back to cartesian system.
% Then doing some data fitting in order to plot the result.

tmp_field(1,:) = linspace(0,max(x_range),floor(x_res/2));
tmp_field(2,:)=zeros(1,length(tmp_field(1,:)));
tmp_field(3,:) = pointSource_field(Phase,T,atomPos,tmp_field(1,:),...
    0,z,lambda,type);
[tmp_field(1,:),tmp_field(2,:)]=cart2pol(tmp_field(1,:),tmp_field(2,:));
for i=0:1:360
    rad_i = deg2rad(i);
    if i>0 && i<360
        field = cat(2,field,[tmp_field(1,2:end)+rad_i;tmp_field(2,2:end);tmp_field(3,2:end)]);
    elseif i==0
        field = tmp_field;
    else
        continue
    end
end

[field(1,:),field(2,:)] = pol2cart(field(1,:),field(2,:));

x_min = min(field(1,:));
x_max = max(field(1,:));
y_min = min(field(2,:));
y_max = max(field(2,:));

% Using a interpolation model in order to plot the field distribution.
F = scatteredInterpolant(field(1,:)',field(2,:)',field(3,:)','natural','none');
[xq,yq] = meshgrid(x_min:((x_max-x_min)/x_res):x_max,...
    y_min:((y_max-y_min)/y_res):y_max);
zq = F(xq,yq);
field=zq;
% figure;
% colormap('jet');
% image(xq(1,:),yq(:,1),real(zq),'CDataMapping','scaled');
% title(['real part at z=',num2str(z),'um']);
% xlabel('x');ylabel('y');
% figure;
% colormap('jet');
% image(xq(1,:),yq(:,1),imag(zq),'CDataMapping','scaled');
% title(['imaginary part at z=',num2str(z),'um']);
% xlabel('x');ylabel('y');
figure;
colormap('jet');
image(xq(1,:),yq(:,1),abs(zq),'CDataMapping','scaled');
title(['absolute value at z=',num2str(z),'um']);
xlabel('x');ylabel('y');
end

