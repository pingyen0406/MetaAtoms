function Phase = NorPhase(inputPhase,normalize)
% This function normalizes the phase data and sets the first one = -1.
% Input:
% Phase: Ascending or descending phase data. Can be either normalized or
% not.
% Output: Normalized phase between [-1,0]. The first one starts from -1 and
% ascends.
arraysize = size(inputPhase);
Phase = zeros(arraysize);
% Normalize the phase 
if normalize==true
    inputPhase=inputPhase/2/pi;
end
% Sbstract the first phase
for i=1:arraysize(1)
    for j=1:arraysize(2)
        if i==1 && j==1
        else
            Phase(i,j)=inputPhase(i,j)-inputPhase(1,1);
        end
        while Phase(i,j)>0
            Phase(i,j)=Phase(i,j)-1;
        end
        while Phase(i,j)<-1
            Phase(i,j)=Phase(i,j)+1;
        end
    end
end
Phase(1,1)=-1;
% Check the value no >1 or <-1
end