
% disp() function with ability to be quiet and suppress output.

function dispq(string, displayswitch)

persistent enabled

if isempty(enabled)
  enabled = true;
end

if nargin == 2
    % either specify with 1 or 0 to show message,
    if isnumeric(displayswitch)
        enabled = displayswitch == 1;
    end
    % or use 'q' or 'quiet':
    if ischar(displayswitch)
       enabled = ~any(strcmpi(displayswitch,{'q','quiet'}));
    end
end

if enabled
  disp(string);
end





