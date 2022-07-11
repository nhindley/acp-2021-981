
% quick function for finding the mean of a load of angles without having to
% deal with wraparound.

% note this isn't perfect, but pretty good. For example, if you put in exactly 
% [0 180] unfortunately it breaks and you can't do it.

%                sum_i_from_1_to_N sin(a[i])
% a = arctangent ---------------------------
%                sum_i_from_1_to_N cos(a[i])


% input in DEGREES, unless additional 'radians' flag is specified.


function OUT = meanangle(varargin)

IN = varargin{1};

if nargin == 2 && strcmpi(varargin{2},'radians')
    OUT = atan2( nansum(sin(IN)) , nansum(cos(IN)) );
else
    OUT = atan2d( nansum(sind(IN)) , nansum(cosd(IN)) );
end
    
















