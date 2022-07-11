

% LETTER FOR AXES.
%
% Input position in NORMALISED UNITS on the axes.
%
%

function t = nph_axesletter(varargin)

ax = gca;

% See if we're doing labels for 3D or 2D axes:
if isnumeric(varargin{1}) && isnumeric(varargin{2}) && isnumeric(varargin{3})
    type = 3;
else
    type = 2;
end


switch type
    case 2
        posx = varargin{1};
        posy = varargin{2};
        posz = 0;
        str = varargin{3};
        varargin = {varargin{4:end}};
    case 3
        posx = varargin{1};
        posy = varargin{2};
        posz = varargin{3};
        str = varargin{4};
        varargin = {varargin{5:end}};
end

xpos = ax.XLim(1) + posx * abs(diff(ax.XLim));
ypos = ax.YLim(1) + posy * abs(diff(ax.YLim));
zpos = ax.ZLim(1) + posz * abs(diff(ax.ZLim));

t = text(gca,xpos,ypos,zpos,str,'HorizontalAlignment','left','VerticalAlignment','bottom',varargin{:});





    
    %     
% switch nargin
%     case 3
% 
% 
% % if you entered a 3D position:
% if isnumeric(str)
%     posz = str;
%     str = varargin{1};
%     varargin{1} = cell(1,1);
%     
%     xpos = ax.XLim(1) + posx * abs(diff(ax.XLim));
%     ypos = ax.YLim(1) + posy * abs(diff(ax.YLim));
%     zpos = ax.YLim(1) + posz * abs(diff(ax.ZLim));
%     
%     t = text(gca,xpos,ypos,zpos,str,'HorizontalAlignment','left','VerticalAlignment','bottom',varargin{:});
%     
% else % normal 2D:
%     
%     xpos = ax.XLim(1) + posx * abs(diff(ax.XLim));
%     ypos = ax.YLim(1) + posy * abs(diff(ax.YLim));
%     
%     t = text(gca,xpos,ypos,str,'HorizontalAlignment','left','VerticalAlignment','bottom',varargin{:});
%     
%     
% end





