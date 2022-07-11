
% draw a line around the current axes:

function borderline(varargin)

axx = gca;

% if ~ishandle(varargin{1})
%     axx = gca;
%     varargin = varargin(2:end);
% end

    hold on; plot(axx.XLim([1 2 2 1 1]),axx.YLim([1 1 2 2 1]),'color','k','linewi',1.5,varargin{:});

end





