


function [Tp,BG] = airs_4dp_detrend(T,Dim,Method)

warning off
  
%reshape the data so only the desired dimension is left
  sz = size(T);
  dims = 1:1:numel(sz);
  order = unique([Dim,dims],'stable');  
  T2 = permute(T,order);
  T2 = reshape(T2,[sz(order(1)),prod(sz(order(2:end)))]);

  %%now we need to remove the 4th order polynomial
  %we can do this with either full polyfit or a massively stripped-down version with no error or quality checking
  % in tests, I see no difference whatsoever, so the fast method is the default (method 2)
  %but retain the option of doing it the old way just in case
  if ~exist('Method','var'); Method = 2; end
  
  %create backround array
  BG = T2.*NaN;  
  
  
  if Method == 1; 
    %"classical" way, with all error checking using built-in

    for iRow=1:1:size(BG,2)
      p = polyval(polyfit(1:1:size(BG,1),T2(:,iRow)',4),1:1:size(BG,1));
      BG(:,iRow) = p;
    end
    
  elseif Method == 2; 
    %fast way - most stuff stripped out. Seems to work...
    %see comment by Matt tearle on
    %https://uk.mathworks.com/matlabcentral/answers/1836-multiple-use-of-polyfit-could-i-get-it-faster
    
    x = repmat((1:1:size(BG,1))',1,prod(sz(order(2:end))));
    n = 4;
    m = size(x,2);
    p = zeros(n+1,m);
    for k = 1:m
      M = repmat(x(:,k),1,n+1);
      M = bsxfun(@power,M,0:n);
      p(:,k) = M\T2(:,k);
    end
    p = flip(p,1);   
    for iRow=1:1:size(BG,2); BG(:,iRow) = polyval(squeeze(p(:,iRow)),1:1:size(BG,1)); end

  end
  
  %back to common code
  
  %put back into original shape
  BG = reshape(BG,[sz(order(1)),(sz(order(2:end)))]);
  BG = permute(BG,order);
 
  %compute perturbations and return
  Tp = T - BG;

return