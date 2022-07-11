

% a nice colorwheel I ripped of the internet somewhere

function cw = colorwheel(varargin)

cw = [
    0 82 161
    0 124 198
    0 156 166
    0 150 60
    135 192 16
    255 227 5
    255 186 18
    255 127 25
    255 77 31
    255 15 33
    231 21 131
    128 35 132] ./ 255;

%%%% allow selecting individual colors:
switch nargin
    case 0
        
    case 1
        if isnumeric(varargin{1})
            cw = cw(varargin{1},:);
        end
end


end

















