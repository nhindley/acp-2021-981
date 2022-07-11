

% Function daynumber.m

% From a matlab datenum, outputs the day of the year.
% OR, specify multiple inputs that could be fed into datenum without
% breaking it

function dayofyear = daynumber(varargin)

switch nargin
    case 1 %%%% JUST MATLAB DATENUM
        datenumber = varargin{1};
        dayofyear = floor(datenumber - datenum(year(datenumber),0,0));
        
    otherwise %%%% YYYY MM DD or similar
        datenumber = datenum(varargin{:});
        dayofyear = floor(datenumber - datenum(year(datenumber),0,0));
        
end










