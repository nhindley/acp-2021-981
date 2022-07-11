
% varargin2struct.m, NPH 20171106
% 
% converts varargin name-pair arguments to a structure where the fields are
% the names and the vars are the field values.
% 
% requires equal number of cell name-pair arguments
% 
% V = varargin2struct(varargin{:});
% 

function V = varargin2struct(varargin)

IN = varargin;

% IN should be a cell object containing varargin inputs

% Parse input:
if isodd(length(IN))    
   error('Incorrect number of name-pair arguments.')
end

% Generate structure:
V = struct;
for i = 1:2:length(IN)/2
    V.(lower(IN{i})) = IN{i+1};
end

% If no output:
if ~nargout
   varargout = V;
end


end



