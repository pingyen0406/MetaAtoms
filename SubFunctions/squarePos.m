function [atomPos_X,atomPos_Y] = squarePos(midPoint,period,size)
% Generating rectangular x,y position of square lattice at given period and 
% size.
% Input:
% midPoint: middle point of the lens, [x,y]
% size: "circle" -> radius; "square" -> [width,height]
% Output: 2*N array

% Genaral concept:
% Create a quarter of the lens, flip it, and then check if it's odd or not.
% If it's odd, we needs to create those points either x or y=0.  
Nx = floor(size(1)/period)+1;
Ny = floor(size(2)/period)+1;
% Create a quarter of the lens
pos_x = -period*(Nx-1)/2:period:period*(Nx-1)/2;
pos_y = -period*(Ny-1)/2:period:period*(Ny-1)/2;
[atomPos_X,atomPos_Y] = meshgrid(pos_x+midPoint(1),pos_y+midPoint(2));


end

