
% 
% nph_convertclass
% 
% Converts the input data [IN] into the class specified by [FORMAT],
% preserving relative values.
% 
% Useful for when you want to convert to a class to save on disk/memory
% space and then, crucially, convert back again to the original values,
% albeit at a slightly reduced precision.
% 
% e.g. convert doubles to arbitrary range of uint8 to save memory, then
% restore to double. You will lose a little bit of precision depending on
% the range of values in [IN].
% 
% INPUTS ==================================================================
%
% IN                    -   data to be converted
%
% OPTION                -   if OPTION is:
% 
%   - a string, it expects it to be a string for the new class, ie convert
%   to this class.
% 
%   - a structure, it expects it to be the Conversion Metadata structure
% created by a previous conversion, ie it will restore to original class
% using this metadata.

% I tried using varargin but I could seem to extract the structure from the
% cell object varargin{1}. Was very weird.

% OUTPUTS:

% OUT                   -   the converted data in the new/original class.
% ConversionMetadata    -   contains data about the conversion so that
%                           it can be restored.


% SUPPORTED CLASSES:
% int8 int16 int32 int64
% uint8 uint16 uint32 uint64

function varargout = nph_convertclass(IN,option,varargin)
    
if exist('varargin','var'), ConvSpec = varargin; end

switch nargout
    
    case 2 % CONVERTING
        % option is a string specifying the desied class
        [OUT,ConvMeta] = nph_convert(IN,option,ConvSpec{:});
        varargout{1} = OUT;
        varargout{2} = ConvMeta;
        
    case 1 % RESTORING
        % Option is a ConvMeta structure containing metadata about the
        % original conversion:
        OUT = nph_restore(IN,option,ConvSpec{:});
        varargout{1} = OUT;
        
end % end switch

end


% Convert:
function [A,ConvMeta] = nph_convert(IN,option,varargin)

classout = char(option);

% Parse Inputs
supportedclasses = {'uint8' 'uint16' 'uint32' 'uint64' ...
    'int8' 'int16' 'int32' 'int64'};

if ~any(strcmpi(classout,supportedclasses))
    disp(['Unsuppoted output class: ' classout])
    return
end

% First, normalise to range of values:
if ~isempty(varargin)
    
    % Check to see if the user has pre-entered the Max and Min (useful if
    % you know what they'll always be in advance so you only need one
    % ConvData structure to convert back many arrays, potentially)
    
    if any(strcmpi('MaxIn',varargin))
        ConvMeta.MaxIn = varargin{find(strcmpi('MaxIn',varargin))+1};
    end
    
    if any(strcmpi('MinIn',varargin))
        ConvMeta.MinIn = varargin{find(strcmpi('MinIn',varargin))+1};
    end
    
else % otherwise do it adhoc
    ConvMeta.MaxIn = max(IN(:));
    ConvMeta.MinIn = min(IN(:));
end


% Subscribe output metadata
ConvMeta.ClassIn = class(IN);
ConvMeta.ClassOut = classout;


% Subtract the minimum to start from zero:
A = (IN - ConvMeta.MinIn);

% Divide by maximum to normalise to one:
A = A ./ (ConvMeta.MaxIn - ConvMeta.MinIn);

% Apply scale for number of values depending on chosen class:
switch classout
    
    % UNSIGNED (much easier!)
    case {'uint8' 'uint16' 'uint32' 'uint64'}
        A = A * abs(double(intmax(classout)));
        eval(['A = ' classout '(A);']);
        
        % Estimate precision out:
        ConvMeta.PrecisionOut = (ConvMeta.MaxIn - ConvMeta.MinIn) / abs(double(intmax(classout)));
        
        % SIGNED
    case {'int8' 'int16' 'int32' 'int64'}
        A = A * (abs(double(intmax(classout))) + abs(double(intmin(classout))));
        A = A - abs(double(intmin(classout)));
        eval(['A = ' classout '(A);']);
        
        % Estimate precision out:
        ConvMeta.PrecisionOut = (ConvMeta.MaxIn - ConvMeta.MinIn) / (abs(double(intmax(classout))) + abs(double(intmin(classout))));
        
    otherwise
        disp(['Unsuppoted output class: ' classout])
        return
        
end



% Sanity check it's worked properly:
ConvMeta.MaxOut = max(A(:));
ConvMeta.MinOut = min(A(:));




end


% Restore:
function A = nph_restore(IN,option,varargin)

ConvMeta = option;

classout = ConvMeta.ClassOut;

% First, convert to double:
A = double(IN);

switch classout
    
    % Normalise to 1:
    case {'uint8' 'uint16' 'uint32' 'uint64'}
        
        A = A / double(intmax(classout));
        
    case {'int8' 'int16' 'int32' 'int64'}
        
        A = A + abs(double(intmin(classout)));
        A = A / double(intmax(['u' classout]));
        
end

% Multiply by the range of original values:
A = A * (ConvMeta.MaxIn - ConvMeta.MinIn);

% And add the min:
A = A + ConvMeta.MinIn;


end

















