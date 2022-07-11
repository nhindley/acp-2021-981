

% Brighten an RGB color triplet or Nx3 colourmap by the factor FACTOR.

% EDIT: Also adjusted to work on images

function OUT = nph_brighten(IN,factor)

% sz = size(IN);
% 
% switch sum(sz ~= 1)
%     case 1
        % single color
        
        OUT = factor * IN;
        
        switch class(IN)
            case 'uint8'
                
            otherwise
                OUT(OUT > 1) = 1;
%                 OUT(OUT > 255) = 255;
                
        end
        
        
%     case 2
%         % Nx3 colormap
%         IN = factor * IN;
%         IN(IN > 1) = 1;
%     case 3
%         % MxNx3 image
%         IN(:,:,2) = factor * IN(:,:,2);
%         
% end

% OUT = hsv2rgb(IN_hsv);

% ensure class is presevered:
OUT = cast(OUT,'like',IN);









