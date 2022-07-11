

% List all user functions required by the input function
% 
% IN is a EITHER
%  - a STRING specifiying the filepath to the function
%  - the STRING name of the function itself
% To specify or bundle up multiple files, enter a cell of strings of
% function names as {'function1.m','function2.m',...}
%
%
% OUT is a cell array of paths to functions 
% 
% OPTIONS:
%  'full' gives the full path to the required functions, otherwise default
%  is simply to list their names.
%  'bundle' bundles all the required functions into a zip folder in the
%  current directory.

%%%% NPH EXTRAS:
% Added an ability to bundle up all the used functions into a zip file.


function OUT = resolveall(IN,varargin)

if isempty(varargin)
    varargin = {''};
end

%%%% single or multiple files...
if ~iscell(IN) % if it's just one, put it in a cell for ease of indexing
    IN = {IN};
end

%%%% now for each file in the cell array...
flist = {};
for i = 1:length(IN)
    
    if ~exist(IN{i},'file') % ie if it's NOT a full file path,
        IN{i} = which(IN{i}); % find the file path
    end
    
    %%%% find all required files and functions...
    disp(['Getting required files and functions for ' IN{i} '...'])
    singleflist = matlab.codetools.requiredFilesAndProducts(IN{i})';
    
    %%%% and cat them up
    flist = cat(1,flist,singleflist);
    
end


%%%% OUTPUT the file list
if any(strcmpi(varargin,'full'))
    OUT = flist;
else
    OUT = cell(size(flist));
    for i = 1:length(flist)
        fun = strsplit(flist{i},'/');
        OUT{i} = fun{end};
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% BUNDLE OPTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% method: copy all used functions to a temporary folder, navigate to it,
% then zip everything in it, then delete the folder.


if any(strcmpi(varargin,'bundle'))
    
    if length(IN) == 1
        outputfilename = IN{1};
    else
        outputfilename = 'multiple_function';
    end
    
    % make temp folder
    tempname = ['temp-' datestr(today)];
    mkdir([pwd '/' tempname])
    % cd([pwd tempname])
    
    % now copy each function to the folder
    for i = 1:length(flist)
%         system(['cp ' flist{i} ' ' pwd '/' tempname]);
        copyfile(flist{i},[pwd '/' tempname]);
    end
    
    % zip it all up:
    system(['zip -qj functionbundle ' pwd '/' tempname '/*']);
    
    % remove the temporary folder
    delete([pwd '/' tempname '/*'])
    rmdir([pwd '/' tempname])
    
    % change name to the function name:
%     system(['mv ' pwd '/functionbundle.zip ' pwd '/' outputfilename '_bundle.zip']);
    movefile([pwd '/functionbundle.zip'],[pwd '/' outputfilename '_bundle.zip']);
    
    
    % display the finished bundle:
%     disp(['All functions bundled into: ' pwd '/functionbundle.zip'])
    disp(['All functions bundled into: ' pwd '/' outputfilename '_functionbundle.zip'])
    
end








