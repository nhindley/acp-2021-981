
% Uses full dispersion relation to calculate a GW's intrinsic frequency from
% its spatial components k, l and m and latitude.

% input WAVELENGTHS (km) please, not wavenumbers. This is so we don't
% confuse between spatial and angular wavenumbers!

% specify Latitude in degrees, please!

function w = gw_intfreq(lx,ly,lz,lat)

% First, Wavenumbers...

k = (2*pi)/lx;

l = (2*pi)/ly;

m = -(2*pi)/lz;

%% Coriolis Parameter...

% convert lat to radians

lat_r = (abs(lat)/360) * (2*pi);

omega = 7.2722*10^-5;

f = 2*omega*sin(lat_r);

H = 7; %km

%% Brunt-Vaisala Frequency...

%N_sq = 0.5*10^-3;    % typical for the stratosphere

N = 0.020943; % rad/s, corresponds to a period of about 5 mins.

%% Dispersion relation from Fritts and Alexander 2003

w = sqrt( ( (N.^2)*(k.^2 + l.^2) + (f.^2)*(m.^2 + 1./(4.*H.^2)) ) ...
        ./ ( k.^2 + l.^2 + m.^2 + 1./(4.*H.^2) ) );

%disp({'rad/s'})


end

% output is RADIANS PER SECOND



