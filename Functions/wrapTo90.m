function lon = wrapTo90(lon)
%wrapTo90 Wrap angle in degrees to [-90 90]
%
%   lonWrapped = wrapTo90(LON) wraps angles in LON, in degrees, to the
%   interval [-90 90] such that 90 maps to 90 and -90 maps to -90.
%   (In general, odd, positive multiples of 90 map to 90 and odd,
%   negative multiples of 90 map to -90.)
%
%   See also wrapTo360, wrapTo2Pi, wrapToPi.

% Copyright 2007-2008 The MathWorks, Inc.

q = (lon < -90) | (90 < lon);
lon(q) = wrapTo360(lon(q) + 90) - 270;