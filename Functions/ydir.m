

function ydir(varargin)

if isempty(varargin)
    set(gca,'ydir','normal');
else
    set(gca,'ydir',varargin{:});
end


end

















