

% Enhanced disp() function that supports multiple colours using Yair
% Altman's cprintf, and also flushes output straight away if we're in
% Octave.




function nph_disp(str,formatspec)

if ~exist('formatspec','var')
    disp(str)
    if isOctave, fflush(stdout); end
    return
end


cprintf(formatspec,[str '\n']);


% If it's Octave, flush to screen now.
% You can also stop Octave buffering screen messages by "more off".
if isOctave, fflush(stdout); end


% cprintf('text',    'regular black text ');
% cprintf('hyper',   'followed %s','by ');
% cprintf('key',     '%d colored ',5);
% cprintf('-comment','& underlined ');
% cprintf('err',     'elements:\n');
% cprintf('cyan',    'cyan ');
% cprintf('_green',  'underlined green ');
% cprintf(-[1,0,1],  'underlined magenta ');
% cprintf([1,0.5,0], 'and multi-\nline orange\n');
% cprintf('*blue',   'and *bold* (R2011b+ only)\n');