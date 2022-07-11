

% For int8, we get +-2^8 (256 total) levels - i.e. -128:127 (including zero)


function doubles = uint82double(integers,max_val,min_val)



doubles = min_val + ((max_val - min_val) * (double(integers) ./ 255));

% doubles = doubles ./ 255;

% doubles = doubles .* (max_val - min_val);

% doubles = doubles + min_val;



% 
% 
% 
% 
% % For int8, we get +-2^8 (256 total) levels - i.e. -128:127 (including zero)
% 
% 
% function [integers] = double2uint8(doubles,max_val,min_val)
% 
% z = doubles - min_val;
% 
% z = z ./ (max_val - min_val);
% 
% z = z .* 255;
% 
% integers = uint8(z);
% % 
% % 
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




% 
% 
% 
% doubles = double(integers);
% 
% doubles = (doubles + 128) ./ 255;
% 
% doubles = doubles .* conversion_factors.norm_factor;
% 
% doubles = doubles + conversion_factors.offset;
% 
% 
% 
% 




















