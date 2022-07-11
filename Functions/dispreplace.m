

% Replaces last line in command window with new string

% Does it by removing last 100 characters of whatever you just displayed.

function  dispreplace(string)

d = '\b';

d = repmat(d,1,length(string));

fprintf(d)

disp(string)









