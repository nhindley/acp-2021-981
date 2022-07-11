

function xdir(varargin)

if isempty(varargin)
    set(gca,'xdir','normal');
else
    set(gca,'xdir',varargin{:});
end


end

















