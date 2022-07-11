
% Originally created by Corwin Wright, University of Bath 2016

% INPUTS
% 
% D1,D2,D3  = X,Y,Z meshgrids on the coordinate frame
% IN        = Data to be plotted
% contour_levels    = array of levels for colour contours
% alpha_levels      = array of levels for alpha transparency
% cmap              = colormap

% OUTPUTS
% 
% OUT       = A 1xlength(contour_levels) cell object containing the patch
%             objects of all the surfaces.
% 

% UPDATE: Now able to output the patch objects as a cell of patch objects,
% one for each isosuface level. This makes it easy to edit properties or
% delete layers etc after generating.


function OUT = nph_isoblobber(D1,D2,D3,IN,contour_levels,alpha_levels,cmap,varargin)

nlevels = length(contour_levels);

if ~exist('cmap','var')
    cmap = flipud(cbrewer('div','RdBu',length(contour_levels)));
end

OUT = cell(1,nlevels);

vis = 'on';

if nargin > 7
    if any(strcmpi(varargin{:},'noplot'))
        vis = 'off';
    end
end

for i = 1:nlevels
    
    if strcmpi(vis,'on'), hold on; end
    
    patchid = patch(isosurface(D1,D2,D3,IN,contour_levels(i)),'visible',vis);
    
    
    try % try the isonormals, if D1,D2,D3 are regular meshgrids:
    if strcmpi(vis,'on'), hold on; end
        isonormals(D1,D2,D3,IN,patchid);
    catch
    end
    
    patchid.FaceColor = cmap(i,:);
    patchid.EdgeColor = 'none';
    patchid.FaceAlpha = alpha_levels(i);
    patchid.BackFaceLighting = 'reverselit';
    
    OUT{i} = patchid;
    
%     drawnow;

end

%% POLISH =================================================================

% set(gcf,'colormap',cmap);
% set(gca,'clim',[contour_levels(1) contour_levels(end)])

lighting gouraud;
material dull;

% camlight;

% if nargin > 7
%     if any(strcmpi(varargin{:},'noplot'))
%         close(gcf);
%     end
% end



