

%%%% don't think this works yet!! :(



%%%% so, this may or may not work. use with caution!
%%%% a WEIGHTED MEAN FOR ANGLES function!


% quick function for finding the mean of a load of angles without having to
% deal with wraparound.

% note this isn't perfect, but pretty good. For example, if you put in exactly 
% [0 180] unfortunately it breaks and you can't do it.

%                sum_i_from_1_to_N sin(a[i] * w[i])
% a = arctangent ---------------------------
%                sum_i_from_1_to_N cos(a[i] * w[i])

% input in DEGREES, unless additional 'radians' flag is specified.

%%%% USAGE:
%%%% OUT = wmeanangle(IN,w)
%%%% OUT = wmeanangle(IN,w,'radians');


function OUT = wmeanangle(varargin)

IN = varargin{1};

if nargin > 1
    if isnumeric(varargin{2})
        w = varargin{2};
    end
end

if any(strcmpi(varargin,'radians'))
    radianflag = 1;
else
    radianflag = 0;
end

if radianflag
    OUT = atan2( nansum(sin( IN .* w)),nansum(cos( IN .* w)));
else
    OUT = atan2d(nansum(sind(IN .* w)),nansum(cosd(IN .* w)));
end
    

