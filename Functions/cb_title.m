



function OUT = cb_title(cbar,str,varargin)

cpos = cbar.Position;
% cpos = cbarposition;

xoffset = 0.01;
yoffset = 0;

w = 4*cpos(3);
h = cpos(4);

x = cpos(1) + (cpos(3)/2) - w/2 + xoffset;
y = cpos(2) + cpos(4) + 0.01 + yoffset;

OUT = annotation('textbox','string',str,varargin{:},...
    'position',[x y w h],'edgecolor','none','fontsize',get(cbar,'fontsize'),...
    'horizontalalignment','center','verticalalignment','bottom');







