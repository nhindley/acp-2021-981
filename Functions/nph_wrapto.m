

% nph_wrapto

% I'm sick and tired of wrapTo90 and stuff not working.

% here's my own version

function OUT = nph_wrapto(IN,val)

val = abs(val);

OUT = mod(IN+val,2.*val)-val;















