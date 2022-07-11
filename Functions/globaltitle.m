

function varargout = globaltitle(varargin)

a = axes('visible','off','units','normalized','position',[0 0 1 0.95]);
t1 = title(varargin{:});
set(t1,'Visible','on');

if nargout
    varargout{1} = t1;
end

end




