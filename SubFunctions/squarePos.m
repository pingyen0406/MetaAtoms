function atomPos = squarePos(type,midPoint,period,size)
% Generating x,y position of square lattice at given period and size.
% Input:
% type: Define shape of the lens. Can be either"circle" or "rect".
% midPoint: middle point of the lens, [x,y]
% size: "circle" -> radius; "square" -> [width,height]
% Output: 2*N array

% Genaral concept:
% Create a quarter of the lens, flip it, and then check if it's odd or not.
% If it's odd, we needs to create those points either x or y=0.  
atomPos = zeros(0,0);

if type=="circle"
    % Create a quarter of the lens
    N = floor(2*size/period)+1;
    for i=1:floor(N/2)    
        for j=1:floor(N/2)
            x_now = period*i;
            y_now = period*j;
            if x_now^2+y_now^2 < size^2
                atomPos=cat(2,atomPos,[x_now;y_now]);
            end               
        end
    end
    % Check it's even or odd and use the corresponding method  
    if mod(N,2)~=0
        atomPos = oddOperation(atomPos,floor(N/2),period); 
    else
        atomPos = evenOperation(atomPos,period);
    end
    
elseif type=="square"
    N = floor(size/period)+1;
    % Create a quarter of the lens
    for i=1:floor(N/2)
        for j=1:floor(N/2)
            x_now = period*i;
            y_now = period*j;
            atomPos=cat(2,atomPos,[x_now;y_now]);     
        end
    end
    % Check it's even or odd and use the corresponding method       
    if mod(N,2)~=0
        atomPos = oddOperation(atomPos,floor(N/2),period);    
    else
        atomPos = evenOperation(atomPos,period);
    end        
else
    error('Wrong input shape');
end
atomPos = atomPos+[midPoint(1);midPoint(2)];
atomPos = sortrows(atomPos')';

end

% For odd numbers, we need to take those either x or y=0 into account.
function newPos = oddOperation(currentPos,N,period)
tmpPos = zeros(2,0);
for i=1:N
    tmpPos(1,i)=period*i;
end
tmpPos = cat(2,tmpPos,-tmpPos);
tmpPos = cat(2,[0;0],tmpPos,flip(tmpPos,1));
newPos = cat(2,tmpPos,currentPos,[-currentPos(1,:);+currentPos(2,:)],...
    [-currentPos(1,:);-currentPos(2,:)],[currentPos(1,:);-currentPos(2,:)]);
end

% For even numbers, just move by -0.5*(period,period) and then copy it.
function newPos = evenOperation(currentPos,period)
currentPos = currentPos-0.5*period;
newPos = cat(2,currentPos,[-currentPos(1,:);+currentPos(2,:)],...
    [-currentPos(1,:);-currentPos(2,:)],[currentPos(1,:);-currentPos(2,:)]);
end