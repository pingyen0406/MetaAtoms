function Phase = NorPhase(Phase)
% This function normalizes the phase data and find the relative phase from 
% the first one.

% Normalize the phase 
if max(max(Phase))>1 | min(min(Phase))<-1
    Phase=Phase/2/pi;
end
% Find the relative value
for i=2:length(Phase)
    Phase(1,i)=Phase(1,i)-Phase(1,1);
end
Phase(1,1)=-1;
% Check the value no >1 or <-1
for i=1:length(Phase)
    while Phase(1,i)>0
        Phase(1,i)=Phase(1,i)-1;
    end
    while Phase(1,i)<-1
        Phase(1,i)=Phase(1,i)+1;
    end
end
end