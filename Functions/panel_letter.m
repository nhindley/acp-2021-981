
% NPH March 2018
% Insert a panel letter into an axes using an annotation textbox.
% Position units are normalized to axes' position in figure.
% Extra inputs are passed to the text command.
% 
% T = panel_letter(axes,posx,posy,str,varargin)
%

function t = panel_letter(ax,posx,posy,str,varargin)

apos = get(ax,'position');

defaultspec = {...
    'units','normalized','fitboxtotext','on', ...
    'edgecolor','none','backgroundcolor',[1 1 1],'facealpha',0, ...
    'horizontalalignment','left','verticalalignment','bottom'};

tpos = [apos(1) + posx*apos(3) apos(2) + posy*apos(4) 0.1 0.1];

t = annotation('textbox','position',tpos,'string',str,defaultspec{:},varargin{:});

end




































