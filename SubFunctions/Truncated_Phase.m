function Phase = Truncated_Phase(Phase,start_index,stop_index) 
% Truncate the phase between the given indices.
if Phase(1,start_index)-Phase(1,stop_index)<0
    Phase = NorPhase(Phase(1,start_index:stop_index));
    Phase(1,end)=Phase(1,end)+1;
else
    Phase = NorPhase(Phase(1,start_index:stop_index));
end
end