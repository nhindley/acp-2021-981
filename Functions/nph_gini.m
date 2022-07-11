
% Compute GINI COEFFICIENT of a vector or N-D array along a specified
% dimension.

% NPH version of Oleg Komarov's code on matlab file exchange.


IN = [0 0 0 0 5 6 7 0 0 0 0];

% function G = nph_gini(IN,varargin)

% Cite Plougonven et al (2013?) for applying this method to GWs and the
% like.

% GINICOEFF Calculate the Gini coefficient, a measure of statistical dispersion 
%   
%   GINICOEFF(IN) Calculate the Gini coeff for each column
%       - In can be a numeric vector/2D matrix of positive real values. 
%   
%   GINICOEFF(...,DIM) Calculate Gini coeff along the dimension selected by DIM
%       - dim: 1 column-wise computation(DEFAULT) 
%              2 row-wise computation
%   
%   GINICOEFF(...,NOSAMPLECORR) Don't apply sample correction
%       - nosamplecorr: true  or 1, apply correction (DEFAULT)
%                       false or 0, don't apply and divide by 'n' 
%
%    [COEFF, IDX] = ...
%       - coeff : n by 1 vector with gini coefficients where n is the number of 
%                 series in IN. The coefficient ranges from 1 (total inequality, 
%                 one unit receives everything) to 0 (total equality, everyone 
%                 receive the same amount).
%       - IDX   : n by 1 logical index. True means that the computation of the gini
%                 coefficient for that series has been skipped due to negative values
%                 or one-element series.
%
% NOTE: NaNs are ignored. The Gini coeff is not computed for those series with 
%       negative values or with just one element. A warning is issued if any series 
%       has negative values or is just a one-element series and IDX is not called 
%       explicitly.
%
% The formula can be expressed as:
% NUM = sum[(n+1-i)*Yi]  for i = 1,...,n; where Yi is monotonically increasing
% DEN = sum[Yi]
% G  = [n+1-2*(NUM/DEN)]*1/n      without sample correction 
% Gc = [n+1-2*(NUM/DEN)]*1/(n-1)  with sample correction (DEFAULT)
%
% Examples:
%   - In = [1,1,1,1,9,3,1,1,1,7,NaN];
%     G  = ginicoeff(In)
%
% Additional features:
% - <a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/26452-gini-coefficient','-browser')">FEX ginicoeff page</a>

% Author: Oleg Komarov (oleg.komarov@hotmail.it) 
% Tested on R14SP3 (7.1) and on R2009b
% 21 jan 2010 - Created 
% 05 feb 2010 - Per Jos (10584) suggestion: sample correction is now the default; if elements in a series < 2 --> NaN. Edited help. Added link to FEX page.
% 15 jun 2010 - NaNs were not ignored due to a misplacing in the lines of code

% CHECK
% ------------------------------------------------------------------------

% NINPUTS
narginchk(1,3);

% IN
if ~isnumeric(IN); error('ginicoeff:fmtIn', 'IN should be numeric'); end

% empty input:
if isempty(IN), coeff = NaN; IDX = NaN; return; end

% DIM
if nargin == 1 || isempty(dim); dim = 1;
elseif ~(isnumeric(dim) && (dim == 1 || dim == 2))
    error('ginicoeff:fmtDim', 'DIM should be 1 or 2');
end

% NOSAMPLECORR
if nargin < 3 || isempty(nosamplecorr)
    nosamplecorr = false;
elseif  ~isscalar(nosamplecorr) || ~(islogical(nosamplecorr) || ismembc(nosamplecorr,[0,1]))
    error('ginicoeff:fmtnosamplecorr','nosamplecorr invalid format');
end


return

%%

IN = rand(1,10);

IN(3) = NaN;

sz = size(IN);
sz = sz(sz ~= 1);
type = length(sz);

% fix the annoying 1xN versus Nx1 problem:
if type == 1
    IN = IN(:);
end
osz = size(IN); % original corrected size in

% assign dimension
% if ~isempty(varargin)
%     dim = varargin{1};
% else
    dim = 1;
% end

ginicoeff(IN)

switch type
    case 0
        return
    otherwise
end


%%

IN = rand(5,10);

% IN = IN(:);

IN(3,3) = NaN;

naninds = isnan(IN);

% no nans:
IN(naninds) = 0;

% no negs:

% Sort In
IN = sort(IN,dim,'ascend');

% Calculate frequencies for each series
freq = flip(cumsum(flip(~naninds,dim),dim),dim);

% Total numel
if dim == 1
    totNum = freq(1,:);
else
    totNum = freq(:,1);
end

% Totals
tot = sum(IN,dim);

% Gini's coefficient
coeff = totNum+1-2*(sum(IN.*freq,dim)./tot);

coeff = coeff./totNum





%%

if all(naninds)
    disp('all nans')
    return
end

nIN = IN;
IN = IN;

naninds = isnan(IN)
IDXnan = isnan(IN)



% use this so as to preserve size:
nIN(naninds) = 0
IN(IDXnan) = 0


nIN = sort(nIN,dim,'ascend')
IN = sort(IN,dim,'ascend')

nfreq = cumsum(~naninds,dim,'reverse')
freq = flip(cumsum(flip(~IDXnan,dim),dim),dim)

ntot = sum(nIN,dim)
tot = sum(IN,dim)

% Total numel
if dim == 1
    totNum = freq(1,:)
    ntotNum = max(nfreq,[],dim)
else
    totNum = freq(:,1)
    ntotNum = max(nfreq,[],dim)
end

coeff = totNum+1-2*(sum(IN.*freq,dim)./tot)
ncoeff = ntotNum+1-2*(sum(nIN.*nfreq,dim)./tot)

coeff = coeff./totNum;
ncoeff = ncoeff./ntotNum;









%%

return


%%

nosamplecorr = 1;

% IN VECTOR
if isvector(IN); IN = IN(:); dim = 1; end 

% ENGINE
% ------------------------------------------------------------------------

% Negative values or one-element series (not admitted)
IDXnan = isnan(IN);
IDX = any(IN < 0,dim) | sum(~IDXnan,dim) < 2;
if dim == 1
    IN(:,IDX) = 0;
else
    IN(IDX,:) = 0;
end
% if nargout ~= 2 && any(IDX)
%     warning('warnginicoeff:negValues','Check IDX for negative values or one-element series')
% end

% Replace NaNs
IN(IDXnan) = 0;

% Sort In
IN = sort(IN,dim,'ascend');

% Calculate frequencies for each series
freq = flip(cumsum(flip(~IDXnan,dim),dim),dim);

% Total numel
if dim == 1
    totNum = freq(1,:);
else
    totNum = freq(:,1);
end

% Totals
tot = sum(IN,dim);

% Gini's coefficient
coeff = totNum+1-2*(sum(IN.*freq,dim)./tot)

% Sample correction
if nosamplecorr
    coeff = coeff./totNum;
else
    coeff = coeff./(totNum-1);
end









