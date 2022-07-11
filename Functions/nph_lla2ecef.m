
% nph edited version of this lla2ecef function.
% this is to differentiate it from the built-in function, which has inputs
% of radians and m etc, which I don't like :)


% LLA2ECEF - convert latitude, longitude, and altitude to
%            earth-centered, earth-fixed (ECEF) cartesian
% 
% USAGE:
% [x,y,z] = lla2ecef(lat,lon,alt)
% 
% IN
% lat = geodetic latitude (deg)
% lon = longitude (deg)
% alt = height above WGS84 ellipsoid (km)
% 
% OUT
% x = ECEF X-coordinate (km)
% y = ECEF Y-coordinate (km)
% z = ECEF Z-coordinate (km)
% 
% Notes: This function assumes the WGS84 model.
%        Latitude is customary geodetic (not geocentric).
% 
% Source: "Department of Defense World Geodetic System 1984"
%         Page 4-4
%         National Imagery and Mapping Agency
%         Last updated June, 2004
%         NIMA TR8350.2
% 
% Michael Kleder, July 2005

function [x,y,z] = nph_lla2ecef(lat,lon,alt)

% WGS84 ellipsoid constants:
a = 6378.137; % km
e = 8.1819190842622e-2;

% intermediate calculation
% (prime vertical radius of curvature)
N = a ./ sqrt(1 - e^2 .* sind(lat).^2);

% results:
x = (N+alt) .* cosd(lat) .* cosd(lon);
y = (N+alt) .* cosd(lat) .* sind(lon);
z = ((1-e^2) .* (N + alt)) .* sind(lat);












