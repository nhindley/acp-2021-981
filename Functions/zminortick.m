



function zminortick(varargin)

switch nargin
    case 1
        axx = gca;
        tix = varargin{1};
    case 2
        axx = varargin{1};
        tix = varargin{2};
end

axx.ZMinorTick = 'on';
axx.ZAxis.MinorTickValues = tix;
