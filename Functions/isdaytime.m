
% a hacked up version of Francois Beauducel's code, IPGP, Paris, France
% for finding sunrise/sunset times, called sunrise.m. It's available on the
% matlab file exchange.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% My version outputs a LOGICAL 1/0 if the sun is up or not for a given
% scalar or vector set of LAT, LON, and TIME (UTC, datenum).

% assume altitude is MSL.

% OUT = isdaytime(LAT,LON,TIME);

% INPUTS: LAT/LON (degrees), TIME (matlab datenum).
% assume altitude is MSL.

% EDIT: I noticed I have to take into account if the sun is up, but on the
% next/previous day.

% % % % %%% EXAMPLE:
% % % % %%% A simple world map of sunshine for current time:
% % % % [LON,LAT] = meshgrid(-180:180,-90:90);
% % % % TIME = now .* ones(size(LAT));
% % % % figure; pcolor(LON,LAT,double(isdaytime(LAT,LON,TIME))); shading flat;


function OUT = isdaytime(varargin)

if any(strcmpi('lars',varargin))
    varargin = varargin(~strcmpi('lars',varargin));
    larsflag = 1;
else
    larsflag = 0;
end
    
switch length(varargin)
    case 3 % lat lon time
        LAT = varargin{1};
        LON = varargin{2};
        TIME = varargin{3};
        ALT = 0; % height above sea level in m.
    case 4 % lat lon time alt
        LAT = varargin{1};
        LON = varargin{2};
        TIME = varargin{3};
        ALT = varargin{4};
end

% Switch to Lars's definition for AIRS 3D:
if larsflag
    OUT = which_airs_retrieval(LON,LAT,TIME);
    return
end

% Starts computation

% number of days since Jan 1st, 2000 12:00 UT
dy = floor(TIME); % this day.
n2000 = dy - datenum(2000,1,1,12,0,0) + 68.184/86400;

% mean solar moon
Js = n2000 - LON./360;

% solar mean anomaly
M = mod(357.5291 + 0.98560028*Js,360);

% center
C = 1.9148*sind(M) + 0.0200*sind(2*M) + 0.0003*sind(3*M);

% ecliptic longitude
lambda = mod(M + C + 180 + 102.9372,360);

% solar transit
Jt = 2451545.5 + Js + 0.0053*sind(M) - 0.0069*sind(2*lambda);

% Sun declination
delta = asind(sind(lambda)*sind(23.44));

% hour angle
omega = acosd( (sind(-0.83 - 2.076*sqrt(ALT)/60) - sind(LAT).*sind(delta))./(cosd(LAT).*cosd(delta)) );

noon = Jt + datenum(2000,1,1,12,0,0) - 2451545;
sset = real(noon + omega/360); % sometimes you get complex, not sure why!
srise = real(noon - omega/360);

thisday = TIME >= srise & TIME <= sset;
nextday = (TIME-1) >= srise & (TIME-1) <= sset;
prevday = (TIME+1) >= srise & (TIME+1) <= sset;
OUT = thisday | nextday | prevday;

end










function Mask = which_airs_retrieval(Lon,Lat,Time)

%identified whether AIRS-3D is using the day or night retrieval, using the same
%logic as Lars' retrieval does itself (shown in C at bottom of function)
%
%Corwin Wright, c.wright@bath.ac.uk, 14/AUG/2018
%
%input is time in Matlab format, latitude, and longitude (array-safe, must all
%be same shape)
%output is a binary mask, with 1 for daytime and 0 for nighttime, of the same
%size as the inputs

%Number of days and fraction with respect to 2000-01-01T12:00Z
D = Time-datenum(2000,1,1,12,0,0);

%Geocentric apparent ecliptic longitude [rad]
g = (357.529 + 0.98560028 .* D) .* pi ./ 180;
q = 280.459 + 0.98564736 .* D;
L = (q + 1.915 .* sin(g) + 0.020 .* sin(2 * g)) .* pi ./ 180;

%Mean obliquity of the ecliptic [rad]
e = (23.439 - 0.00000036 * D) .* pi ./ 180;

% Declination [rad]
dec = asin(sin(e) .* sin(L));

% Right ascension [rad]...
ra = atan2(cos(e) .* sin(L), cos(L));
 
% Greenwich Mean Sidereal Time [h]...
GMST = 18.697374558 + 24.06570982441908 .* D;

%Local Sidereal Time [h]... 
LST = GMST + Lon ./ 15;

% Hour angle [rad]... 
h = LST ./ 12 .* pi - ra;
 
% Convert latitude... 
Lat = Lat .* pi ./ 180;

%hence, solar zenith angle [deg].
SZAd = acos(sin(Lat) .* sin(dec) + cos(Lat) .* cos(dec) .* cos(h)) .* 180 / pi;

%create mask
Mask = SZAd .*NaN;
Mask(SZAd <  96) = 1; %day
Mask(SZAd >= 96) = 0; %night



end


