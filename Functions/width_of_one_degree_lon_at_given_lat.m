

% INPUT LATITUDE IN DEGREES, output is the width of one degree longitude at
% this latitude in km, assuming spherical earth.

function w = width_of_one_degree_lon_at_given_lat(lat)

% fix the explosion at the poles...

lat(lat == 90)  =  89.999;

lat(lat == -90) = -89.999;

% put lat in radians:

latr = (lat/360)*(2*pi);

Re = 6371; % km, approx radius of the earth.

w = ((2*pi) / 360) * (Re*cos(latr));

% output is in km.

end




















