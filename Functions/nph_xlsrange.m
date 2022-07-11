

% I've only been working on it a few hours but I'm already sick of having
% to specificy excel spreadsheet rows and columns using numbrs AND letters
% when exporting to them.

% This function takes NUMERIC row and column indeces and converts it to a
% string range such as 'A1:E5' etc.


function OUT = nph_xlsrange(startcol,startrow,endcol,endrow)

OUT = [upper(num2str(alphabet(startcol))) num2str(startrow) ':' upper(num2str(alphabet(endcol))) num2str(endrow)];








