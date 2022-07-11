

% NPH 20171005 UoB
% 
% List all user functions required by the input function
% 
% IN is a EITHER
%  - a STRING specifiying the filepath to the function
%  - the STRING name of the function itself
% 
% OUT is a cell array of paths to functions 
% 
% OPTIONS:
%  'full' gives the full path to the required functions, otherwise default
%  is simply to list their names.
% 



function OUT = requiredfunctions(IN,varargin)

if isempty(varargin)
    varargin = {''};
end

if ~exist(IN,'file') % ie if it's NOT a full file path,
    IN = which(IN); % find the file path
end

flist = matlab.codetools.requiredFilesAndProducts(IN)';


if any(strcmpi(varargin{:},'full'))
    
    OUT = flist;
    
else
    
    OUT = cell(size(flist));
    
    for i = 1:length(flist)
        fun = strsplit(flist{i},'/');
        OUT{i} = fun{end};
    end
    
end












