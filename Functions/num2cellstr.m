

% Convert numeric vector to a cell array of strings of the numbers.
% Very useful for setting axes/colorbar labels from vectors.

% N-D inputs are linearised.

function OUT = num2cellstr(IN)

OUT = cellstr(num2str(IN(:)))';







