
% move_object(handle,direction,distance)

function move_object(handle,varargin)

dirs = varargin(1:2:end-1);
dists = varargin(2:2:end);

for i = 1:length(dirs)
    
    direction = dirs{i};
    distance = dists{i};

pos = get(handle,'position');

if strcmp(direction,'up')
    
    pos(2) = pos(2) + distance;
    
end

if strcmp(direction,'down')
    
    pos(2) = pos(2) - distance;
    
end

if strcmp(direction,'left')
    
    pos(1) = pos(1) - distance;
    
end

if strcmp(direction,'right')
    
    pos(1) = pos(1) + distance;
    
end


set(handle,'position',pos);

end



