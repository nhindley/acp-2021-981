

% For int8, we get +-2^8 (256 total) levels - i.e. -128:127 (including zero)


function [integers] = double2uint8(doubles,max_val,min_val)

% integers = single(doubles) - min_val; % If you're going to uint8, you won'd miss double precision ;)

% integers = integers - min_val;

integers = 255 * ((single(doubles) - min_val) ./ (max_val - min_val));

% integers = z .* 255;

% integers = uint8(z);
% 
% 
% a = doubles - min(doubles(:));
% conversion_factors.offset = min(doubles(:));
% 
% 
% a = a ./ max(a(:));
% conversion_factors.norm_factor = max(a(:));
% 
% a = (255 .* a) - 128;
% 
% 
% if max(a(:)) > 127 || min(a(:)) < -128,
%     
%     disp('Problem with converting to int8.')
%     
% end
% 
% integers = int8(a);



% if abs(min(doubles(:))) <= abs(max(doubles(:))),
%     
%     scale_factor = max(doubles(:));
%     
% else
%     
%     scale_factor = min(doubles(:));
%     
% end

% integers = int8( 127 * (doubles ./ scale_factor) );

























