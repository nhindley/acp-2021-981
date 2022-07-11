
% quickly calculate coriolis parameter (inertial frequency) 
% at any given latitude...

function f = inertialfreq(lat)


%% Coriolis Parameter...

% convert lat to radians

lat_r = (abs(lat)/360) * (2*pi);

omega = 7.2722*10^-5; % rotational frequency of the earth (rad/s)

f = 2*omega*sin(lat_r); % rad/s

end


