
% quickly calculate coriolis parameter (inertial frequency) 
% at any given latitude...

function t = inertialperiod(lat)


%% Coriolis Parameter...

% convert lat to radians

lat_r = (abs(lat)/360) * (2*pi);

omega = 7.2722*10^-5; % rotational frequency of the earth (rad/s)

f = 2*omega*sin(lat_r); % rad/s

t = (1./f)*(2*pi); % in seconds....

t = (t/3600); % in HOURS

%disp({'Hours'})

end


