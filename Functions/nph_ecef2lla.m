

% nph edited version of this ecef2lla function.
% this is to differentiate it from the built-in function, which has inputs
% of radians and m etc, which I don't like :)



% ECEF2LLA - convert earth-centered earth-fixed (ECEF)
%            cartesian coordinates to latitude, longitude,
%            and altitude
%
% USAGE:
% [lat,lon,alt] = ecef2lla(x,y,z)
%
% IN
% x = ECEF X-coordinate (km)
% y = ECEF Y-coordinate (km)
% z = ECEF Z-coordinate (km)
% 
% OUT
% lat = geodetic latitude (deg)
% lon = longitude (deg)
% alt = height above WGS84 ellipsoid (km)
% 
%
% Notes: (1) This function assumes the WGS84 model.
%        (2) Latitude is customary geodetic (not geocentric).
%        (3) Inputs may be scalars, vectors, or matrices of the same
%            size and shape. Outputs will have that same size and shape.
%        (4) Tested but no warranty; use at your own risk.
%        (5) Michael Kleder, April 2006

function [lat,lon,alt] = nph_ecef2lla(x,y,z)

% % Put it in m:
% x = x .* 1000; y = y .* 1000; z = z .* 1000; 

% WGS84 ellipsoid constants:
a = 6378.137;
e = 8.1819190842622e-2;

% calculations:
b   = sqrt(a^2*(1-e^2));
ep  = sqrt((a^2-b^2)/b^2);
p   = sqrt(x.^2+y.^2);
th  = atan2d(a*z,b*p);
lon = atan2d(y,x);
lat = atan2d((z+ep^2.*b.*sind(th).^3),(p-e^2.*a.*cosd(th).^3));
N   = a./sqrt(1-e^2.*sind(lat).^2);
alt = p./cosd(lat)-N;

% return lon in range [0,2*pi)
% lon = mod(lon,2*pi);
lon = wrapTo180(lon);

% correct for numerical instability in altitude near exact poles:
% (after this correction, error is about 2 millimeters, which is about
% the same as the numerical precision of the overall function)

k=abs(x)<1 & abs(y)<1;
alt(k) = abs(z(k))-b;

% lat = r2d(lat);
% lon = r2d(lon);
% alt = alt ./ 1000;






