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
      
elseif symmetry==true  
    field = Focal_Rotate(Phase,T,atomPos,x_range,y_range,z,x_res,y_res,lambda);  
     
else
    error("Wrong input on symmetry");
end

end
