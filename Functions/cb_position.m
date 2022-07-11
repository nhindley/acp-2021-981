

% get nice colorbar positon

function OUT = cb_position(axes_handle,horz_gap)

if nargin == 1
    horz_gap = 0.05;
end

apos = get(axes_handle,'Position');

x = apos(1) + apos(3) + horz_gap;
y = apos(2) + 0.25*apos(4);

w = 0.025;
h = 0.5*apos(4);

OUT = [x y w h];










