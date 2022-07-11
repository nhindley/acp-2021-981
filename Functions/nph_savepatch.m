

% If you want to save a patch object and plot it later using patch(), you
% need to remove the type, parent, annotation and being deleted fields.
% 
% P = nph_savepatch(handle_to_patch_object);
% figure; patch(P);
% 
% or:
% nph_savepatch('savepath',handle_to_patch_object);
% 

function varargout = nph_savepatch(varargin)

switch nargin
    case 1
        IN = varargin{1};
    case 2
        savepath = varargin{1};
        IN = varargin{2};
end

P = struct;

f = fieldnames(IN);
            
notfields = {'type','parent','annotation','beingdeleted'};

for i = 1:length(f)
    
    if ~any(strcmpi(f{i},notfields))
        
        P.(f{i}) = IN.(f{i});
        
    end
    
end

switch nargin
    case 1
        varargout{1} = P;
    case 2
        save(savepath,'P')
end
        

% str = str(1:end-1); % removre trailing comma

% %% change to lat lon alt:
% 
% % SGWEX MODEL CENTRE POINT:
% clat = -54.5269874182436;
% clon = -37.1367495437688;
% 
% [lat,lon] = reckon(clat,clon,km2deg(quadadd(S.XData,S.YData)),atan2d(S.XData,S.YData));
% 
% 

% 
% %% PLOT ON NEW, CORRECT LAT LON GRID 
% 
% S.XData = lon;
% S.YData = lat;

% eval(['figure; hold on; patch(' str ')']);

% grid on;











