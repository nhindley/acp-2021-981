

% FIXED the below font/line scaling bug with the function sdf() which is
% copied at the end of the function below :)

% Another bug fix: recently, fontsizes and scaling seem to be not working
% correctly. You can fix this under Export Settings in the figure menu,
% where I've set the font scaleing to be 100% of the original value - i.e.
% none. This should have been saved now, but I can't seem to do it
% programmatically so do be aware of the problem.

% Quick fix repalcement for export_fig, which doesn't work on
% Yosemite/2014b due to Ghostscript not working properly. A known problem
% but doesn't seem to be solution yet.

% Exporting as a PDF via saveas was fine, but you always ended up with an
% A4 page with big white borders.

% This approach sets the physical size of the paper (in inches) to be what
% ever fraction of the screen is occupied by the figure. If the figure is
% fullscreen, then the paper size is set to the physical width in inces of
% my dell screen. If half the size of the screen, the paper is that size.

% -> get size of screen
% -> scale pixels to inches using scale_factor or default
% -> set paper size to be size of figure paper
% -> set paper position to be full size of figure paper.
% -> save using saveas.

% Could also implement a border trimming aspect by having an input trim
% vector [left_trim right_trim top_trim bottom_trim] with values from 0 to
% 1, where we set a fractional amount to trim off the end of the
% paperpositon parameter. But that can wait :)

% Another important note is that the fontsize and linewidth properties are
% tied up with the size of the page. This means you can't arbitrarily screw
% around with the scale_factor without screwing up your fonts, I believe.
% At least I think that's one of the problems I encountered.



%%

function nph_export_fig(filename,type,scale_factor)

% useable screen size on dell monitor = 8 X 20.5 inches / 1200 X 1920 pix
% useable screen size on macbook monitor = 7.75 X 14.5 inches.

% It's best to set the scale factor to be of order the size of your
% monitor, ideally the width, since the width is what I'm using to do the
% pixels to inches conversion. This ensures that the font sizes etc that
% you see on the screen will be exactly that size when saved.

% Set default scaling factor...
if nargin <= 2, scale_factor = 25; end;


% figpos = get(gcf,'position');
% screensize = get(0,'screensize');
% 
% px2in = scale_factor/screensize(3);
% 
% figheight = figpos(4)*px2in;
% figwidth  = figpos(3)*px2in;
% 
% set(gcf,'papersize',[figwidth figheight])
% set(gcf,'paperposition',[0 0 figwidth figheight])


% v THIS SEEMS TO BE CRITICAL NOW FOR FONTS AND LINES SCALING
% sdf(gcf,'default');
sdf(gcf);


saveas(gcf,filename,type)






%% WE NOW USE SDF BELOW ===================================================

function sdf(varargin)
% SDF Set the line width and fonts of a figure
% 
% sdf(fig)
% 
% where fig is the figure number. If the figure number is omitted, the 
% currently active figure is updated. Edit the file to set you own style 
% settings.
%
% sdf(fig, 'stylename')
% applies a pre-configured style from the File-->Export Setup menu of the
% figure's window. The stylename should be one of the 'Export Styles'
% section of the dialog.
%
% The function allows applying the same settings as through the 
% File-->Export Setup-->Apply menu of the figure, but much faster and 
% without the annoying clicking. 
%
% Example
%   figure(1);      t=0:0.1:10;   plot(t, sin(t));
%   sdf(1)
%   pause
%   sdf(1,'PowerPoint')

% Andrey Popov, Hamburg, 2009

%% Parse the input data
default = true;
if nargin == 0       % no input - take current fig and apply default style
    fig = gcf();
else                 % at least 1 input
    if ischar(varargin{1})  % style name
        default = false;
        style_name = varargin{1};
        fig = gcf();
    else
        fig = varargin{1};
        figure(fig);        % just in case it does not exist
        if nargin == 2
            default = false;
            style_name = varargin{2};
        end
    end
end

%% Apply a style
if default      % Apply a default style
    style = struct();
    style.Version = '1';
    style.Format = 'eps';
    style.Preview = 'none';
    style.Width = 'auto';
    style.Height = 'auto';
    style.Units = 'centimeters';
    style.Color = 'rgb';
    style.Background = 'w';          % '' = no change; 'w' = white background
%     style.FixedFontSize = '20';
%     style.ScaledFontSize = 'auto';
%     style.FontMode = 'fixed';
%     style.FontSizeMin = '8'; % commented out as seemed to rescale upon execution
%     style.FixedLineWidth = '1';
%     style.ScaledLineWidth = 'auto';
%     style.LineMode = 'fixed';
%     style.LineWidthMin = '0.5';
    style.FontName = 'auto';
    style.FontWeight = 'auto';
    style.FontAngle = 'auto';
    style.FontEncoding = 'latin1';
    style.PSLevel = '2';
    style.Renderer = 'auto';
    style.Resolution = 'auto';
    style.LineStyleMap = 'none';
    style.ApplyStyle = '0';
    style.Bounds = 'loose';
    style.LockAxes = 'on';
    style.ShowUI = 'on';
    style.SeparateText = 'off';

    hgexport(fig,'temp_dummy',style,'applystyle', true);

else    % Apply an existing style, defined as in the Export dialog
    % The files are in folder   fullfile(prefdir(0),'ExportSetup');
    style = hgexport('readstyle',style_name);
    hgexport(fig,'temp_dummy',style,'applystyle', true);
end























